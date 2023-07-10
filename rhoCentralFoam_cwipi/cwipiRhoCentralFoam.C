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

Application
    rhoCentralFoam

Group
    grpCompressibleSolvers

Description
    Density-based compressible flow solver based on
    central-upwind schemes of Kurganov and Tadmor with
    support for mesh-motion and topology changes.

\*---------------------------------------------------------------------------*/

#include "cwipiPstream.H"
#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote(
        "Density-based compressible flow solver based on"
        " central-upwind schemes of Kurganov and Tadmor with"
        " support for mesh-motion and topology changes.");

#define NO_CONTROL
#include "postProcess.H"
#include "addCheckCaseOptions.H"
#include "setRootCaseLists.H"
#include "createTime.H"
#include "createDynamicFvMesh.H"
#include "createFields.H"
#include "createFieldRefs.H"
#include "createTimeControls.H"

    turbulence->validate();

#include "cwipiBaseFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "readFluxScheme.H"

    const dimensionedScalar v_zero(dimVolume / dimTime, Zero);

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info
        << "\nStarting time loop\n"
        << endl;

    if (cwipiSwitch == true)
    {
#include "cwipiInitialise.H"
        while (runTime.run())
        {
#include "cwipiUpdateFields.H"
#include "cwipiRunTime.H"
#include "rhoCentralFoam.H"
            runTime.write();
#include "cwipiTimeStep.H"
            runTime.printExecutionTime(Info);
        }
    }
    else
    {
        while (runTime.run())
        {
#include "cwipiUpdateFields.H"
#include "rhoCentralFoam.H"
            runTime.write();
            runTime.printExecutionTime(Info);
        }
    }

    Info << "End\n"
         << endl;

    return 0;
}

// ************************************************************************* //
