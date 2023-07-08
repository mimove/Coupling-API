#include "readTimeControls.H"

if (!LTS)
{
#include "setDeltaT.H"

    ++runTime;

    // Do any mesh changes
    mesh.update();
}

// --- Directed interpolation of primitive fields onto faces
surfaceScalarField rho_pos(interpolate(rho, pos));
surfaceScalarField rho_neg(interpolate(rho, neg));
surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));
volScalarField rPsi("rPsi", 1.0 / psi);
surfaceScalarField rPsi_pos(interpolate(rPsi, pos, T.name()));
surfaceScalarField rPsi_neg(interpolate(rPsi, neg, T.name()));
surfaceScalarField e_pos(interpolate(e, pos, T.name()));
surfaceScalarField e_neg(interpolate(e, neg, T.name()));
surfaceVectorField U_pos("U_pos", rhoU_pos / rho_pos);
surfaceVectorField U_neg("U_neg", rhoU_neg / rho_neg);
surfaceScalarField p_pos("p_pos", rho_pos *rPsi_pos);
surfaceScalarField p_neg("p_neg", rho_neg *rPsi_neg);
surfaceScalarField phiv_pos("phiv_pos", U_pos &mesh.Sf());

// Note: extracted out the orientation so becomes unoriented
phiv_pos.setOriented(false);
surfaceScalarField phiv_neg("phiv_neg", U_neg &mesh.Sf());
phiv_neg.setOriented(false);

// Make fluxes relative to mesh-motion
if (mesh.moving())
{
    surfaceScalarField meshPhi(mesh.phi());
    meshPhi.setOriented(false);
    phiv_pos -= meshPhi;
    phiv_neg -= meshPhi;
}

volScalarField c("c", sqrt(thermo.Cp() / thermo.Cv() * rPsi));
surfaceScalarField cSf_pos(
    "cSf_pos",
    interpolate(c, pos, T.name()) * mesh.magSf());
surfaceScalarField cSf_neg(
    "cSf_neg",
    interpolate(c, neg, T.name()) * mesh.magSf());
surfaceScalarField ap(
    "ap",
    max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero));
surfaceScalarField am(
    "am",
    min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero));
surfaceScalarField a_pos("a_pos", ap / (ap - am));
surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));
surfaceScalarField aSf("aSf", am *a_pos);

if (fluxScheme == "Tadmor")
{
    aSf = -0.5 * amaxSf;
    a_pos = 0.5;
}

surfaceScalarField a_neg("a_neg", 1.0 - a_pos);

phiv_pos *= a_pos;
phiv_neg *= a_neg;

surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);

// Reuse amaxSf for the maximum positive and negative fluxes estimated by the central scheme
amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

if (mesh.nInternalFaces())
{
    scalarField sumAmaxSf(fvc::surfaceSum(amaxSf)().primitiveField());

    CoNum = 0.5 * gMax(sumAmaxSf / mesh.V().field()) * runTime.deltaTValue();

    meanCoNum = 0.5 * (gSum(sumAmaxSf) / gSum(mesh.V().field())) * runTime.deltaTValue();
}

Info << "Mean and max Courant Numbers = "
     << meanCoNum << " " << CoNum << endl;

if (LTS)
{
    {
        volScalarField &rDeltaT = trDeltaT.ref();

        scalar rDeltaTSmoothingCoeff(
            runTime.controlDict().getOrDefault<scalar>(
                "rDeltaTSmoothingCoeff",
                0.02));

        // Set the reciprocal time-step from the local Courant number
        rDeltaT.ref() = max(
            1 / dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
            fvc::surfaceSum(amaxSf)()() / ((2 * maxCo) * mesh.V()));

        // Update tho boundary values of the reciprocal time-step
        rDeltaT.correctBoundaryConditions();

        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

        Info << "Flow time scale min/max = "
             << gMin(1 / rDeltaT.primitiveField())
             << ", " << gMax(1 / rDeltaT.primitiveField()) << endl;
    }

    ++runTime;
}

Info << "Time = " << runTime.timeName() << nl << endl;

phi = aphiv_pos * rho_pos + aphiv_neg * rho_neg;

surfaceVectorField phiU(aphiv_pos *rhoU_pos + aphiv_neg * rhoU_neg);
// Note: reassembled orientation from the pos and neg parts so becomes oriented
phiU.setOriented(true);

surfaceVectorField phiUp(phiU + (a_pos * p_pos + a_neg * p_neg) * mesh.Sf());

surfaceScalarField phiEp(
    "phiEp",
    aphiv_pos *(rho_pos *(e_pos + 0.5 * magSqr(U_pos)) + p_pos) + aphiv_neg * (rho_neg * (e_neg + 0.5 * magSqr(U_neg)) + p_neg) + aSf * p_pos - aSf * p_neg);

// Make flux for pressure-work absolute
if (mesh.moving())
{
    surfaceScalarField meshPhi(mesh.phi());
    meshPhi.setOriented(false);
    phiEp += meshPhi * (a_pos * p_pos + a_neg * p_neg);
}

volScalarField muEff("muEff", turbulence->muEff());
volTensorField tauMC("tauMC", muEff *dev2(Foam::T(fvc::grad(U))));

// --- Solve density
solve(fvm::ddt(rho) + fvc::div(phi));

// --- Solve momentum
solve(fvm::ddt(rhoU) + fvc::div(phiUp));

U.ref() = rhoU() / rho();
U.correctBoundaryConditions();
rhoU.boundaryFieldRef() == rho.boundaryField() * U.boundaryField();

if (!inviscid)
{
    solve(fvm::ddt(rho, U) - fvc::ddt(rho, U) - fvm::laplacian(muEff, U) - fvc::div(tauMC));
    rhoU = rho * U;
}

// --- Solve energy
surfaceScalarField sigmaDotU(
    "sigmaDotU",
    (
        fvc::interpolate(muEff) * mesh.magSf() * fvc::snGrad(U) + fvc::dotInterpolate(mesh.Sf(), tauMC)) &
        (a_pos * U_pos + a_neg * U_neg));

solve(
    fvm::ddt(rhoE) + fvc::div(phiEp) - fvc::div(sigmaDotU));

e = rhoE / rho - 0.5 * magSqr(U);
e.correctBoundaryConditions();
thermo.correct();
rhoE.boundaryFieldRef() == rho.boundaryField() * (e.boundaryField() + 0.5 * magSqr(U.boundaryField()));

if (!inviscid)
{
    solve(fvm::ddt(rho, e) - fvc::ddt(rho, e) - fvm::laplacian(turbulence->alphaEff(), e));
    thermo.correct();
    rhoE = rho * (e + 0.5 * magSqr(U));
}

p.ref() = rho() / psi();
p.correctBoundaryConditions();
rho.boundaryFieldRef() == psi.boundaryField() * p.boundaryField();

// Blend fields to reference solution
U = (U * acousticBlending) - (U_inf * (acousticBlending - 1));
thermo.rho() = (thermo.rho() * acousticBlending) - (rho_inf * (acousticBlending - 1));
thermo.p() = (thermo.p() * acousticBlending) - (p_inf * (acousticBlending - 1));
thermo.T() = (thermo.T() * acousticBlending) - (T_inf * (acousticBlending - 1));

turbulence->correct();
