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

#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "volFields.H"
#include "psiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    class cwipi
    {
    public:
        // Default constructor
        cwipi(
            const Foam::Time &runTime,
            const fvMesh &mesh,
            const psiThermo &thermo);

        // Default destructor
        ~cwipi();

        // Create coupling
        void createCoupling(
            const fvMesh &mesh);

        // Update source fields to latest time step
        inline void updateSources(
            const volPointInterpolation &pInterp,
            const psiThermo &thermo)
        {
            LPrime_ = L() - LMean();
            TPrime_ = thermo.T() - TMean();
            sPrime_ = s() - sMean();
            ThetaPrime_ = (TPrime() * fvc::grad(sMean())) - (sPrime() * fvc::grad(TMean()));
            F_u_ = (entropyGradientSwitch * ThetaPrime()) - (lambVectorSwitch * LPrime()) * sourceDamping();
            F_p_ = entropyDerivativeSwitch * ((rhoMean() / thermo.Cp()) * (fvc::ddt(sPrime()) + (UMean() & fvc::grad(sPrime())))) * sourceDamping();
        }

        // Send fields
        void send(
            const pointScalarField &F_0_p,
            const pointVectorField &F_0_u,
            const fvMesh &mesh);

        // Update Lamb vector, entropy and speed of sound
        inline void updateFields(
            const volVectorField &U,
            const psiThermo &thermo)
        {
            L_ = fvc::curl(U) ^ U;
            s_ = (thermo.Cp() * log(thermo.T() / Foam::constant::standard::Tstd)) - ((Foam::constant::physicoChemical::R / thermo.W()) * log(thermo.p() / Foam::constant::standard::Pstd));
            c_ = Foam::sqrt(thermo.gamma() * thermo.p() / thermo.rho());
        }

        // Advance time step
        void updateTime();

        // Access to cwipi switch
        bool isActive()
        {
            return isActive_;
        };

        // Access to continuity sources
        inline volScalarField F_p()
        {
            return F_p_;
        }

        // Access to momentum sources
        inline volVectorField F_u()
        {
            return F_u_;
        }

        // Evaluate send condition
        inline bool sendNow()
        {
            if (cwipiTimeStep == cwipiStep)
            {
                return true;
            }
            else
            {
                return false;
            }
        }

    private:
        // Cwipi switch
        const bool isActive_;

        // Send tag and status
        int sendTag;
        int status;

        // Spatial dimension
        const uint8_t cwipiDim;

        // Parameters for coupling
        const double lambVectorSwitch;
        const double entropyGradientSwitch;
        const double entropyDerivativeSwitch;
        const Foam::label cwipiStep;
        Foam::label cwipiTimeStep;

        // Requisite fields
        const volVectorField UMean_;
        const volScalarField rhoMean_;
        const volVectorField LMean_;
        const volScalarField sMean_;
        const volScalarField cMean_;
        const volScalarField TMean_;
        const volScalarField sourceDamping_;

        volVectorField L_;
        volScalarField s_;
        volScalarField c_;

        volVectorField LPrime_;
        volVectorField ThetaPrime_;
        volScalarField sPrime_;
        volScalarField TPrime_;

        volScalarField F_p_;
        volVectorField F_u_;

        // Point fields
        // pointScalarField F_0_p_;
        // pointVectorField F_0_u_;

        // List of argument names
        const char *sourceFieldNames;

        // Mesh info for cwipi
        std::vector<double> pointCoords = std::vector<double>(1);
        std::vector<int> connecIdx = std::vector<int>(1);
        std::vector<int> connec = std::vector<int>(1);

        // Cwipi source array
        std::vector<double> fieldsToSend = std::vector<double>(1);

        // Vol point interpolation
        // volPointInterpolation pInterp;

        // Access to vol fields
        inline const volVectorField UMean()
        {
            return UMean_;
        }

        inline const volScalarField rhoMean()
        {
            return rhoMean_;
        }

        inline const volVectorField LMean()
        {
            return LMean_;
        }

        inline const volScalarField sMean()
        {
            return sMean_;
        }

        inline const volScalarField cMean()
        {
            return cMean_;
        }

        inline const volScalarField TMean()
        {
            return TMean_;
        }

        inline const volScalarField sourceDamping()
        {
            return sourceDamping_;
        }

        inline volVectorField L()
        {
            return L_;
        };

        inline volScalarField s()
        {
            return s_;
        };

        inline volScalarField c()
        {
            return c_;
        };

        inline volVectorField LPrime()
        {
            return LPrime_;
        }

        inline volVectorField ThetaPrime()
        {
            return ThetaPrime_;
        }

        inline volScalarField sPrime()
        {
            return sPrime_;
        }

        inline volScalarField TPrime()
        {
            return TPrime_;
        }
    };

    class fvMesh;

}

// ************************************************************************* //
