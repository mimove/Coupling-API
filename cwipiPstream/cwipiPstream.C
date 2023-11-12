/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "cwipiPstream.H"
#include <cwipi.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // Constructor
    cwipiPstream::cwipiPstream(
        const Foam::Time &runTime,
        const fvMesh &mesh,
        const psiThermo &thermo,
        const cwipiFields &sourceFields,
        const volPointInterpolation &pInterp)
        : sendTag(0),                                                                                               // Set send tag to 0
          status(0),                                                                                                // Set status to 0
          cwipiDim(readInt(runTime.controlDict().lookup("cwipiDim"))),                                              // Get dimension
          lambVectorSwitch(static_cast<Foam::scalar>(readBool(runTime.controlDict().lookup("cwipiLambVector")))),   // Cast lamb vector coefficient to scalar
          entropyGradientSwitch(static_cast<Foam::scalar>(readBool(runTime.controlDict().lookup("cwipiEntropy")))), // Cast entropy gradient coefficient to scalar
          entropyDerivativeSwitch(static_cast<Foam::scalar>(readBool(runTime.controlDict().lookup("cwipiDsDt")))),  // Cast entropy material derivative coefficient to scalar
          cwipiStep(readInt(runTime.controlDict().lookup("cwipiStep"))),                                            // Get time step
          cwipiTimeStep(readInt(runTime.controlDict().lookup("cwipiStep"))),                                        // Assign time step
          UMean_(IOobject("UMean", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh),
          rhoMean_(IOobject("rhoMean", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh),
          LMean_(IOobject("LMean", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh),
          sMean_(IOobject("sMean", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh),
          cMean_(IOobject("cMean", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh),
          TMean_(IOobject("TMean", runTime.timeName(), mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh),
          sourceDamping_(IOobject("sourceDamping", runTime.timeName(), mesh, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE), mesh, dimensionedScalar("one", dimless, 1)),
          LPrime_(IOobject("LPrime", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionSet(0, 1, -2, 0, 0, 0, 0)),
          ThetaPrime_(IOobject("ThetaPrime", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionSet(0, 1, -2, 0, 0, 0, 0)),
          sPrime_(IOobject("sPrime", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionSet(0, 2, -2, -1, 0, 0, 0)),
          TPrime_(IOobject("TPrime", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE), mesh, dimensionSet(0, 0, 0, 1, 0, 0, 0)),
          F_p_(IOobject("F_p_", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh, dimensionSet(1, -3, -1, 0, 0, 0, 0)),
          F_u_(IOobject("F_u_", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh, dimensionSet(0, 1, -2, 0, 0, 0, 0)),
          mesh_(mesh),
          thermo_(thermo),
          sourceFields_(sourceFields),
          pInterp_(pInterp),
          F_0_p_(pInterp_.interpolate(F_p_)),
          F_0_u_(pInterp_.interpolate(F_u_))
    {
        // Catch incorrect parameter at start, no need to throw in updateSources()
        // Also assign sending field names
        switch (cwipiDim)
        {
        case 2:
            sourceFieldNames = "F_0_p,F_0_u,F_0_v";
            Info << "Coupling enabled with 2 physical dimensions" << endl;
            break;
        case 3:
            sourceFieldNames = "F_0_p,F_0_u,F_0_v,F_0_w";
            Info << "Coupling enabled with 3 physical dimensions" << endl;
            break;
        default:
            throw std::invalid_argument("Variable cwipiDim should be 2 or 3.");
            break;
        }

        // Resize mesh vectors to fit
        pointCoords.resize(3 * mesh_.nPoints());
        connecIdx.resize(mesh_.nCells() + 1);
        connec.resize(mesh_.nCells() * 8);
        fieldsToSend.resize((cwipiDim + 1) * mesh_.nPoints());

        // Create mesh connectivity list
        forAll(mesh_.points(), i)
        {
            pointCoords[3 * i + 0] = mesh_.points()[i].x();
            pointCoords[3 * i + 1] = mesh_.points()[i].y();
            pointCoords[3 * i + 2] = mesh_.points()[i].z();
        }
        connecIdx[0] = 0;
        forAll(mesh_.cells(), i)
        {
            connecIdx[i + 1] = connecIdx[i] + 8;
        }
        forAll(mesh_.cells(), i)
        {
            forAll(mesh_.cellShapes()[i], j)
            {
                connec[8 * i + j] = mesh_.cellShapes()[i][j] + 1;
            }
        }

        // Add local control parameters
        cwipi_add_local_int_control_parameter(
            "nSendVars",
            cwipiDim + 1);
        cwipi_add_local_string_control_parameter(
            "sendFieldNames",
            sourceFieldNames);
        cwipi_add_local_int_control_parameter(
            "receiveTag",
            sendTag);

        // Create coupling
        cwipi_create_coupling(
            "cwipiFoam",
            CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
            "FOAM_APE",
            3,
            1,
            CWIPI_STATIC_MESH,
            CWIPI_SOLVER_CELL_VERTEX,
            0,
            "EnSight Gold",
            "text");

        // Sync and dump
        cwipi_synchronize_control_parameter(
            "FOAM_APE");
        cwipi_dump_application_properties();

        // Get send tag
        sendTag = cwipi_get_distant_int_control_parameter(
            "FOAM_APE",
            "receiveTag");

        // Define mesh and locate interpolation
        cwipi_define_mesh(
            "cwipiFoam",
            mesh_.nPoints(),
            mesh_.nCells(),
            pointCoords.data(),
            connecIdx.data(),
            connec.data());
        cwipi_locate(
            "cwipiFoam");
    }

    // Destructor
    cwipiPstream::~cwipiPstream()
    {
    }

    // Send fields
    void cwipiPstream::send()
    {
        // Update source fields
        updateSources();

        // Send data
        cwipi_issend(
            "cwipiFoam",
            "ex1",
            sendTag,
            cwipiDim + 1,
            1,
            0,
            sourceFieldNames,
            fieldsToSend.data(),
            &status);

        // Handle exchange status
        switch (status)
        {
        case CWIPI_EXCHANGE_OK:
            Info << "Exchange Ok" << endl;
            break;
        case CWIPI_EXCHANGE_BAD_RECEIVING:
            Info << "Bad receiving" << endl;
            break;
        default:
            Info << "Error: bad exchange status" << endl;
            break;
        }
    }

    // Advance time step
    void cwipiPstream::updateTime()
    {
        // Compare time step to step
        if (cwipiTimeStep == cwipiStep)
        {
            // Wait
            cwipi_wait_issend(
                "cwipiFoam",
                status);
        }

        // Use modulo operator (think this is correct?)
        cwipiTimeStep = cwipiTimeStep % cwipiStep;

        // Advance time step
        cwipiTimeStep = cwipiTimeStep + 1;
    }

}

// ************************************************************************* //
