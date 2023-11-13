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

#include "fvCFD.H"
#include "volPointInterpolation.H"
#include "volFields.H"
#include "psiThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    class cwipiSwitch
    {
    public:
        // Default constructor
        cwipiSwitch(
            const Foam::Time &runTime)
            : isActive_(readBool(runTime.controlDict().lookup("cwipiSwitch")))
        {
            if (isActive() == true)
            {
                Info << "CWIPI coupling enabled" << endl;
            }
            else
            {
                Info << "CWIPI coupling disabled" << endl;
            }
        };

        // Default destructor
        ~cwipiSwitch(){};

        // Return access to switch
        inline Foam::Switch isActive()
        {
            return isActive_;
        }

    private:
        // Switch
        const Foam::Switch isActive_;
    };

    class cwipiFields
    {
    public:
        // Default constructor
        cwipiFields(
            const fvMesh &mesh,
            const Foam::Time &runTime,
            const volVectorField &U,
            const psiThermo &thermo)
            : L_(IOobject("L", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh, dimensionSet(0, 1, -2, 0, 0, 0, 0)),
              s_(IOobject("s", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh, dimensionSet(0, 2, -2, -1, 0, 0, 0)),
              c_(IOobject("c", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE), mesh, dimensionSet(0, 1, -1, 0, 0, 0, 0)),
              U_(U),
              thermo_(thermo)
        {
            Info << "Cwipi source fields allocated" << endl;
        };

        // Default destructor
        ~cwipiFields(){};

        // Update Lamb vector, entropy and speed of sound
        inline void update()
        {
            // Update Lamb vector
            L_ = fvc::curl(U()) ^ U();

            // Update entropy
            s_ = (thermo_.Cp() * log(thermo_.T() / Foam::constant::standard::Tstd)) - ((Foam::constant::physicoChemical::R / thermo_.W()) * log(thermo_.p() / Foam::constant::standard::Pstd));

            // Update speed of sound
            c_ = Foam::sqrt(thermo_.gamma() * thermo_.p() / thermo_.rho());
        }

        // Access to Lamb vector by reference
        inline volVectorField const &L() const
        {
            return L_;
        }

        // Access to entropy by reference
        inline volScalarField const &s() const
        {
            return s_;
        }

        // Access to speed of sound by reference
        inline volScalarField const &c() const
        {
            return c_;
        }

    private:
        // Lamb vector
        volVectorField L_;

        // Entropy
        volScalarField s_;

        // Speed of sound
        volScalarField c_;

        // Const reference to velocity
        const volVectorField &U_;

        // Const reference to thermo
        const psiThermo &thermo_;

        // Velocity
        inline const volVectorField U() const
        {
            return U_;
        }
    };

    class cwipiPstream
    {
    public:
        // Default constructor
        cwipiPstream(
            const Foam::Time &runTime,
            const fvMesh &mesh,
            const psiThermo &thermo,
            const cwipiFields &sourceFields,
            const volPointInterpolation &pInterp);

        // Default destructor
        ~cwipiPstream();

        // Send fields
        void send();

        // Advance time step
        void updateTime();

        // Access to continuity sources
        inline const volScalarField &F_p() const
        {
            return F_p_;
        }

        // Access to momentum sources
        inline const volVectorField &F_u() const
        {
            return F_u_;
        }

        // Evaluate send condition
        inline const Foam::Switch sendNow() const
        {
            return (cwipiTimeStep == cwipiStep);
        }

    private:
        // Update source fields to latest time step
        inline void updateSources()
        {
            // Compute Lamb vector
            LPrime_ = sourceFields_.L() - LMean();

            // Compute entropy gradient
            TPrime_ = thermo_.T() - TMean();
            sPrime_ = sourceFields_.s() - sMean();
            ThetaPrime_ = (TPrime() * fvc::grad(sMean())) - (sPrime() * fvc::grad(TMean()));

            // Compute continuity equation source
            F_p_ = entropyDerivativeSwitch * ((rhoMean() / thermo_.Cp()) * (fvc::ddt(sPrime()) + (UMean() & fvc::grad(sPrime())))) * sourceDamping();

            // Compute momentum equation source
            F_u_ = (entropyGradientSwitch * ThetaPrime()) - (lambVectorSwitch * LPrime()) * sourceDamping();

            // Pointwise interpolation of sources
            F_0_p_ = pInterp_.interpolate(F_p());
            F_0_u_ = pInterp_.interpolate(F_u());

            // Assign fieldsToSend for either 2D or 3D case
            if (dim() == 2)
            {
                forAll(mesh_.points(), i)
                {
                    fieldsToSend[((dim() + 1) * i) + 0] = F_0_p()[i];
                    fieldsToSend[((dim() + 1) * i) + 1] = F_0_u()[i].x();
                    fieldsToSend[((dim() + 1) * i) + 2] = F_0_u()[i].y();
                }
            }
            else
            {
                forAll(mesh_.points(), i)
                {
                    fieldsToSend[((dim() + 1) * i) + 0] = F_0_p()[i];
                    fieldsToSend[((dim() + 1) * i) + 1] = F_0_u()[i].x();
                    fieldsToSend[((dim() + 1) * i) + 2] = F_0_u()[i].y();
                    fieldsToSend[((dim() + 1) * i) + 3] = F_0_u()[i].z();
                }
            }
        }

        // Send tag and status
        int sendTag;
        int status;

        // Spatial dimension
        const uint8_t cwipiDim;

        // Parameters for coupling
        const scalar lambVectorSwitch;
        const scalar entropyGradientSwitch;
        const scalar entropyDerivativeSwitch;
        const label cwipiStep;
        label cwipiTimeStep;

        // Requisite fields
        const volVectorField UMean_;
        const volScalarField rhoMean_;
        const volVectorField LMean_;
        const volScalarField sMean_;
        const volScalarField cMean_;
        const volScalarField TMean_;
        const volScalarField sourceDamping_;
        volVectorField LPrime_;
        volVectorField ThetaPrime_;
        volScalarField sPrime_;
        volScalarField TPrime_;
        volScalarField F_p_;
        volVectorField F_u_;

        // Mesh
        const fvMesh &mesh_;

        // Thermo
        const psiThermo &thermo_;

        // Cwipi fields object
        const cwipiFields &sourceFields_;

        // Pointwise interpolation
        const volPointInterpolation &pInterp_;

        // Point fields
        pointScalarField F_0_p_;
        pointVectorField F_0_u_;

        // List of argument names
        const char *sourceFieldNames;

        // Mesh info for cwipi
        std::vector<scalar> pointCoords = std::vector<scalar>();
        std::vector<int> connecIdx = std::vector<int>();
        std::vector<int> connec = std::vector<int>();

        // Cwipi source array
        std::vector<scalar> fieldsToSend = std::vector<scalar>();

        // Physical dimension
        inline uint8_t dim() const
        {
            return cwipiDim;
        }

        // Access to vol fields
        inline const volVectorField &UMean() const
        {
            return UMean_;
        }

        inline const volScalarField &rhoMean() const
        {
            return rhoMean_;
        }

        inline const volVectorField &LMean() const
        {
            return LMean_;
        }

        inline const volScalarField &sMean() const
        {
            return sMean_;
        }

        inline const volScalarField &cMean() const
        {
            return cMean_;
        }

        inline const volScalarField &TMean() const
        {
            return TMean_;
        }

        inline const volScalarField &sourceDamping() const
        {
            return sourceDamping_;
        }

        inline const volVectorField &LPrime() const
        {
            return LPrime_;
        }

        inline const volVectorField &ThetaPrime() const
        {
            return ThetaPrime_;
        }

        inline const volScalarField &sPrime() const
        {
            return sPrime_;
        }

        inline const volScalarField &TPrime() const
        {
            return TPrime_;
        }

        inline const pointScalarField &F_0_p() const
        {
            return F_0_p_;
        }

        inline const pointVectorField &F_0_u() const
        {
            return F_0_u_;
        }
    };

}

// ************************************************************************* //
