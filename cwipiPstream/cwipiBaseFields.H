// Requisite CWIPI fields
volVectorField L(
    IOobject(
        "L",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    fvc::curl(U) ^ U);
volScalarField s(
    IOobject(
        "s",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    (thermo.Cp() * log(thermo.T() / Foam::constant::standard::Tstd)) - ((Foam::constant::physicoChemical::R / thermo.W()) * log(thermo.p() / Foam::constant::standard::Pstd)));
volScalarField c(
    IOobject(
        "c",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    Foam::sqrt(thermo.gamma() * thermo.p() / thermo.rho()));

// Get switch for cwipi
const bool cwipiSwitch = readBool(runTime.controlDict().lookup("cwipiSwitch"));
