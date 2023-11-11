/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
     \\/     M anipulation  |
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
    cwipi::cwipi(
        const Foam::Time &runTime,
        const fvMesh &mesh,
        const psiThermo &thermo)
        : isActive_(readBool(runTime.controlDict().lookup("cwipiSwitch"))),                                         // Get switch for cwipi
          sendTag(0),                                                                                               // Set send tag to 0
          status(0),                                                                                                // Set status to 0
          cwipiDim(readInt(runTime.controlDict().lookup("cwipiDim"))),                                              // Get dimension
          lambVectorSwitch(static_cast<Foam::scalar>(readBool(runTime.controlDict().lookup("cwipiLambVector")))),   // Cast lamb vector coefficient to scalar
          entropyGradientSwitch(static_cast<Foam::scalar>(readBool(runTime.controlDict().lookup("cwipiEntropy")))), // Cast entropy gradient coefficient to scalar
          entropyDerivativeSwitch(static_cast<Foam::scalar>(readBool(runTime.controlDict().lookup("cwipiDsDt")))),  // Cast entropy material derivative coefficient to scalar
          cwipiStep(readInt(runTime.controlDict().lookup("cwipiStep"))),                                            // Get time step
          cwipiTimeStep(readInt(runTime.controlDict().lookup("cwipiStep"))),                                        // Assign time step
          UMean_(
              IOobject(
                  "UMean",
                  runTime.timeName(),
                  mesh,
                  IOobject::MUST_READ,
                  IOobject::AUTO_WRITE),
              mesh),
          rhoMean_(
              IOobject(
                  "rhoMean",
                  runTime.timeName(),
                  mesh,
                  IOobject::MUST_READ,
                  IOobject::AUTO_WRITE),
              mesh),
          LMean_(
              IOobject(
                  "LMean",
                  runTime.timeName(),
                  mesh,
                  IOobject::MUST_READ,
                  IOobject::AUTO_WRITE),
              mesh),
          sMean_(
              IOobject(
                  "sMean",
                  runTime.timeName(),
                  mesh,
                  IOobject::MUST_READ,
                  IOobject::AUTO_WRITE),
              mesh),
          cMean_(
              IOobject(
                  "cMean",
                  runTime.timeName(),
                  mesh,
                  IOobject::MUST_READ,
                  IOobject::AUTO_WRITE),
              mesh),
          TMean_(
              IOobject(
                  "TMean",
                  runTime.timeName(),
                  mesh,
                  IOobject::MUST_READ,
                  IOobject::AUTO_WRITE),
              mesh),
          sourceDamping_(
              IOobject(
                  "sourceDamping",
                  runTime.timeName(),
                  mesh,
                  IOobject::READ_IF_PRESENT,
                  IOobject::AUTO_WRITE),
              mesh,
              dimensionedScalar(
                  "one",
                  dimless,
                  1)),
          L_(
              IOobject(
                  "L",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE),
              mesh,
              dimensionSet(0, 1, -2, 0, 0, 0, 0)),
          s_(
              IOobject(
                  "s",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE),
              mesh,
              dimensionSet(0, 2, -2, -1, 0, 0, 0)),
          c_(
              IOobject(
                  "c",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::AUTO_WRITE),
              mesh,
              dimensionSet(0, 1, -1, 0, 0, 0, 0)),
          LPrime_(
              IOobject(
                  "LPrime",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE),
              mesh,
              dimensionSet(0, 1, -2, 0, 0, 0, 0)),
          ThetaPrime_(
              IOobject(
                  "ThetaPrime",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE),
              mesh,
              dimensionSet(0, 1, -2, 0, 0, 0, 0)),
          sPrime_(
              IOobject(
                  "sPrime",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE),
              mesh,
              dimensionSet(0, 2, -2, -1, 0, 0, 0)),
          TPrime_(
              IOobject(
                  "TPrime",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE),
              mesh,
              dimensionSet(0, 0, 0, 1, 0, 0, 0)),
          F_p_(
              IOobject(
                  "F_p_",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE),
              mesh,
              dimensionSet(1, -3, -1, 0, 0, 0, 0)),
          F_u_(
              IOobject(
                  "F_u_",
                  runTime.timeName(),
                  mesh,
                  IOobject::NO_READ,
                  IOobject::NO_WRITE),
              mesh,
              dimensionSet(0, 1, -2, 0, 0, 0, 0))
    {
        if (isActive() == true)
        {
            switch (cwipiDim)
            {
            case 2:
                Info << "Coupling enabled with 2 physical dimensions" << endl;
                break;
            case 3:
                Info << "Coupling enabled with 3 physical dimensions" << endl;
                break;
            default:
                throw std::invalid_argument("Variable cwipiDim should be 2 or 3.");
                break;
            }
        }
    }

    // Destructor
    cwipi::~cwipi()
    {
    }

    // Establish coupling
    void cwipi::createCoupling(const fvMesh &mesh)
    {
        // Resize mesh vectors to fit
        pointCoords.resize(3 * mesh.nPoints());
        connecIdx.resize(mesh.nCells() + 1);
        connec.resize(mesh.nCells() * 8);
        fieldsToSend.resize((cwipiDim + 1) * mesh.nPoints());

        // Create mesh connectivity list
        forAll(mesh.points(), i)
        {
            pointCoords[3 * i + 0] = mesh.points()[i].x();
            pointCoords[3 * i + 1] = mesh.points()[i].y();
            pointCoords[3 * i + 2] = mesh.points()[i].z();
        }
        connecIdx[0] = 0;
        forAll(mesh.cells(), i)
        {
            connecIdx[i + 1] = connecIdx[i] + 8;
        }
        forAll(mesh.cells(), i)
        {
            forAll(mesh.cellShapes()[i], j)
            {
                connec[8 * i + j] = mesh.cellShapes()[i][j] + 1;
            }
        }

        // Assign field names for either 2 or 3d case, otherwise error
        switch (cwipiDim)
        {
        case 2:
            sourceFieldNames = "F_0_p,F_0_u,F_0_v";
            break;
        case 3:
            sourceFieldNames = "F_0_p,F_0_u,F_0_v,F_0_w";
            break;
        default:
            throw std::invalid_argument("Variable cwipiDim should be 2 or 3.");
            break;
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
            mesh.nPoints(),
            mesh.nCells(),
            pointCoords.data(),
            connecIdx.data(),
            connec.data());
        cwipi_locate(
            "cwipiFoam");
    }

    // Send fields
    void cwipi::send(
        const pointScalarField &F_0_p,
        const pointVectorField &F_0_u,
        const fvMesh &mesh)
    {
        // Assign fieldsToSend for either 2 or 3d case, otherwise error
        switch (cwipiDim)
        {
        case 2:
            forAll(mesh.points(), i)
            {
                fieldsToSend[((cwipiDim + 1) * i) + 0] = F_0_p()[i];
                fieldsToSend[((cwipiDim + 1) * i) + 1] = F_0_u()[i].x();
                fieldsToSend[((cwipiDim + 1) * i) + 2] = F_0_u()[i].y();
            }
            break;
        case 3:
            forAll(mesh.points(), i)
            {
                fieldsToSend[((cwipiDim + 1) * i) + 0] = F_0_p()[i];
                fieldsToSend[((cwipiDim + 1) * i) + 1] = F_0_u()[i].x();
                fieldsToSend[((cwipiDim + 1) * i) + 2] = F_0_u()[i].y();
                fieldsToSend[((cwipiDim + 1) * i) + 3] = F_0_u()[i].z();
            }
            break;
        default:
            throw std::invalid_argument("Variable cwipiDim should be 2 or 3.");
            break;
        }

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

    void cwipi::updateTime()
    {
        // Compare time step to step
        if (cwipiTimeStep == cwipiStep)
        {
            // Wait
            Info << "Before wait." << endl;
            cwipi_wait_issend("cwipiFoam", status);
            Info << "After wait." << endl;

            // Reset counter to 0
            cwipiTimeStep = 0;
        }

        // Advance time step
        cwipiTimeStep = cwipiTimeStep + 1;
    }

}

// ************************************************************************* //
