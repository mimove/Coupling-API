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

void cwipiCoupling(const fvMesh& mesh)
{
    scalar nCells = mesh.nCells();
    scalar nPoints = mesh.nPoints();

    double* pointCoords = new double[3*mesh.nPoints()];
    forAll(mesh.points(),i)
    {
        pointCoords[3*i+0]=mesh.points()[i].x();
        pointCoords[3*i+1]=mesh.points()[i].y();
        pointCoords[3*i+2]=mesh.points()[i].z();
    }

    int* connecIdx = new int[mesh.nCells()+1];
    connecIdx[0]=0;
    forAll(mesh.cells(),i)
    {
        connecIdx[i+1]=connecIdx[i]+8;
    }


    int* connec = new int[mesh.nCells()*8];
    forAll(mesh.cells(),i)
    {


        forAll(mesh.cellShapes()[i],j)
        {
            connec[8*i+j]=mesh.cellShapes()[i][j]+1;
        }
    }

    cwipi_add_local_int_control_parameter("nSendVars",3);
    cwipi_add_local_string_control_parameter("sendFieldNames", "u0,v0,w0");
    cwipi_create_coupling("cwipiFoam",
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                          "FOAM_APE",
                          3,
                          1.0,
                          CWIPI_STATIC_MESH,
                          CWIPI_SOLVER_CELL_VERTEX,
                          1,
                          "Ensight Gold",
                          "text");
    cwipi_synchronize_control_parameter("FOAM_APE");
    cwipi_dump_application_properties();
    sendTag = cwipi_get_distant_int_control_parameter("FOAM_APE", "receiveTag");
    cwipi_define_mesh("cwipiFoam", nPoints, nCells, pointCoords,connecIdx,connec);
    cwipi_locate("cwipiFoam");

}


void cwipiSend(const fvMesh& mesh, const volVectorField& vf)
{

    volPointInterpolation pInterp(mesh);
    pointVectorField pL = pInterp.interpolate(vf);

    double* fieldsToSend = new double[3*mesh.nPoints()];
    forAll(mesh.points(),i)
    {
        fieldsToSend[3*i+0]=pL[i].x();
        fieldsToSend[3*i+1]=pL[i].y();
        fieldsToSend[3*i+2]=pL[i].z();
    }


    cwipi_issend("cwipiFoam","ex1",sendTag,3,1,0,"u0,v0,w0",fieldsToSend,&status);

    switch(status)
    {
        case CWIPI_EXCHANGE_OK :
        printf("Exchange Ok\n");
        break;
        case CWIPI_EXCHANGE_BAD_RECEIVING :
        printf("Bad receiving\n");
        break;
        default :
        printf("Error : bad exchange status\n");
    }
}

void cwipiWait()
{
    cwipi_wait_issend("cwipiFoam",status);
}

}

// ************************************************************************* //
