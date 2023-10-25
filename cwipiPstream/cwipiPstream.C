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

    static int sendTag;
    static int status;

    const char *cwipiCoupling(
        const fvMesh &mesh,
        const uint8_t cwipiDim)
    {

        std::string cwipiArgumentList_;

        // Have to allocate these on the heap (stack too small)
        // double *pointCoords = new double[3 * mesh.nPoints()];
        // int *connecIdx = new int[mesh.nCells() + 1];
        // int *connec = new int[mesh.nCells() * 8];

        double pointCoords[3 * mesh.nPoints()];
        int pointCoords[mesh.nCells() + 1];
        int pointCoords[mesh.nCells() * 8];

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

        switch (cwipiDim)
        {
        case 2:
            cwipiArgumentList_ = "F_0_p,F_0_u,F_0_v";
            break;
        case 3:
            cwipiArgumentList_ = "F_0_p,F_0_u,F_0_v,F_0_w";
            break;
        default:
            throw std::invalid_argument("Variable cwipiDim should be 2 or 3.");
            break;
        }
        const char *cwipiArgumentList = cwipiArgumentList_.c_str();
        cwipi_add_local_int_control_parameter("nSendVars", cwipiDim + 1);
        cwipi_add_local_string_control_parameter("sendFieldNames", cwipiArgumentList);
        cwipi_add_local_int_control_parameter("receiveTag", sendTag);

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

        cwipi_synchronize_control_parameter("FOAM_APE");
        cwipi_dump_application_properties();
        sendTag = cwipi_get_distant_int_control_parameter("FOAM_APE", "receiveTag");
        cwipi_define_mesh("cwipiFoam", mesh.nPoints(), mesh.nCells(), pointCoords, connecIdx, connec);
        cwipi_locate("cwipiFoam");

        // delete[] pointCoords;
        // delete[] connecIdx;
        // delete[] connec;

        return cwipiArgumentList;
    }

    void cwipiSend(
        const double *sourceArray,
        const char *cwipiArgumentList,
        const uint8_t cwipiDim)
    {
        cwipi_issend(
            "cwipiFoam",
            "ex1",
            sendTag,
            cwipiDim + 1,
            1,
            0,
            cwipiArgumentList,
            sourceArray,
            &status);

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

    void cwipiWait()
    {
        cwipi_wait_issend("cwipiFoam", status);
    }

    void cwipiReshapeSourceArrays(
        const fvMesh &mesh,
        double *sourceArray,
        const pointScalarField &acousticContinuityEquation,
        const pointVectorField &acousticMomentumEquation,
        const uint8_t cwipiDim)
    {
        switch (cwipiDim)
        {
        case 2:
            forAll(mesh.points(), i)
            {
                sourceArray[((cwipiDim + 1) * i) + 0] = acousticContinuityEquation[i];
                sourceArray[((cwipiDim + 1) * i) + 1] = acousticMomentumEquation[i].x();
                sourceArray[((cwipiDim + 1) * i) + 2] = acousticMomentumEquation[i].y();
            }
            break;
        case 3:
            forAll(mesh.points(), i)
            {
                sourceArray[((cwipiDim + 1) * i) + 0] = acousticContinuityEquation[i];
                sourceArray[((cwipiDim + 1) * i) + 1] = acousticMomentumEquation[i].x();
                sourceArray[((cwipiDim + 1) * i) + 2] = acousticMomentumEquation[i].y();
                sourceArray[((cwipiDim + 1) * i) + 3] = acousticMomentumEquation[i].z();
            }
            break;
        default:
            Info << "Variable cwipiDim should be 2 or 3." << endl;
            break;
        }
    }

}

// ************************************************************************* //
