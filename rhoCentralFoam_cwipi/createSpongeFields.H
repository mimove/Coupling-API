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
