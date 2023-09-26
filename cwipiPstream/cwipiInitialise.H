// Momentum terms
const volVectorField UMean(
    IOobject(
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);
const volScalarField rhoMean(
    IOobject(
        "rhoMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);
const volScalarField pMean(
    IOobject(
        "pMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

// Lamb vector terms
const volVectorField LMean(
    IOobject(
        "LMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);
volVectorField LPrime(
    IOobject(
        "LPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    L - LMean);

// Entropy terms
const volScalarField sMean(
    IOobject(
        "sMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);
volScalarField sPrime(
    IOobject(
        "sPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    s - sMean);

// Temperature terms
const volScalarField TMean(
    IOobject(
        "TMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);
volScalarField TPrime(
    IOobject(
        "TPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    T - TMean);

// Speed of sound terms
const volScalarField cMean(
    IOobject(
        "cMean",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

// Entropy gradient term
volVectorField ThetaPrime(
    IOobject(
        "ThetaPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    (TPrime * fvc::grad(sMean)) - (sPrime * fvc::grad(TMean)));

// Source damping coefficient
const volScalarField sourceDamping(
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
        1));

const bool lambVectorSwitch = readBool(runTime.controlDict().lookup("cwipiLambVector"));
const bool entropyGradientSwitch = readBool(runTime.controlDict().lookup("cwipiEntropy"));
// const bool entropyDerivativeSwitch = readBool(runTime.controlDict().lookup("cwipiDsDt"));
const uint8_t cwipiDim = readInt(runTime.controlDict().lookup("cwipiDim"));
const Foam::label cwipiStep = readInt(runTime.controlDict().lookup("cwipiStep"));
Foam::label cwipiTimeStep = cwipiStep;

volPointInterpolation pInterp(mesh);

volVectorField momentumEquation(
    IOobject(
        "momentumEquation",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    fvc::curl(U) ^ U);
// volScalarField continuityEquation(
//     IOobject(
//         "continuityEquation",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE),
//     (thermo.rho() / thermo.Cp()) * (fvc::ddt(s) + (U & fvc::grad(s))));

pointVectorField acousticMomentumEquation = pInterp.interpolate(momentumEquation);
// pointScalarField acousticContinuityEquation = pInterp.interpolate(continuityEquation);
forAll(mesh.points(), i)
{
    acousticMomentumEquation[i].x() = 0;
    acousticMomentumEquation[i].y() = 0;
    acousticMomentumEquation[i].z() = 0;
    // acousticContinuityEquation[i] = 0;
}

double sourceArray[cwipiDim * mesh.nPoints()];
// double sourceArray[mesh.nPoints()];
const char *cwipiArgumentList = cwipiCoupling(mesh, cwipiDim);
