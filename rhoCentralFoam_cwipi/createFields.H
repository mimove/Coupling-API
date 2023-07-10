// Constants for blending
const dimensionedScalar T_inf(
    "T_inf",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    readScalar(runTime.controlDict().lookup("T_inf")));
const dimensionedScalar p_inf(
    "p_inf",
    dimensionSet(1, -1, -2, 0, 0, 0, 0),
    readScalar(runTime.controlDict().lookup("p_inf")));
const dimensionedScalar rho_inf(
    "rho_inf",
    dimensionSet(1, -3, 0, 0, 0, 0, 0),
    readScalar(runTime.controlDict().lookup("rho_inf")));
const dimensionedVector U_inf(
    "U_inf",
    dimensionSet(0, 1, -1, 0, 0, 0, 0),
    vector(runTime.controlDict().lookup("U_inf")));

// Blending coefficient
const volScalarField acousticBlending(
    IOobject(
        "acousticBlending",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar(
        "one",
        dimless,
        1));

#include "createRDeltaT.H"

Info << "Reading thermophysical properties\n"
     << endl;

autoPtr<psiThermo> pThermo(
    psiThermo::New(mesh));
psiThermo &thermo = pThermo();

volScalarField &e = thermo.he();

Info << "Reading field U\n"
     << endl;
volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

U = (U * acousticBlending) - (U_inf * (acousticBlending - 1));
thermo.rho() = (thermo.rho() * acousticBlending) - (rho_inf * (acousticBlending - 1));
thermo.p() = (thermo.p() * acousticBlending) - (p_inf * (acousticBlending - 1));
thermo.T() = (thermo.T() * acousticBlending) - (T_inf * (acousticBlending - 1));

volScalarField rho(
    IOobject(
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    thermo.rho());

volVectorField rhoU(
    IOobject(
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    thermo.rho() * U);

volScalarField rhoE(
    IOobject(
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    thermo.rho() * (e + 0.5 * magSqr(U)));

surfaceScalarField pos(
    IOobject(
        "pos",
        runTime.timeName(),
        mesh),
    mesh,
    dimensionedScalar("pos", dimless, 1.0));

surfaceScalarField neg(
    IOobject(
        "neg",
        runTime.timeName(),
        mesh),
    mesh,
    dimensionedScalar("neg", dimless, -1.0));

surfaceScalarField phi("phi", fvc::flux(rhoU));

Info << "Creating turbulence model\n"
     << endl;
autoPtr<compressible::turbulenceModel> turbulence(
    compressible::turbulenceModel::New(
        rho,
        U,
        phi,
        thermo));
