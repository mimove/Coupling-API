// Momentum terms
volVectorField UMean(
    IOobject(
        "UMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    U);
volScalarField rhoMean(
    IOobject(
        "rhoMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    thermo.rho());

// Lamb vector terms
volVectorField LMean(
    IOobject(
        "LMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    fvc::curl(U) ^ U);

// Entropy terms
volScalarField sMean(
    IOobject(
        "sMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    s);
volScalarField sPrime(
    IOobject(
        "sPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    s);

// Temperature terms
volScalarField TMean(
    IOobject(
        "TMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    thermo.T());
volScalarField TPrime(
    IOobject(
        "TPrime",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    thermo.T());

// Speed of sound terms
volScalarField cMean(
    IOobject(
        "cMean",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    Foam::sqrt(thermo.gamma() * thermo.p() / thermo.rho()));

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
const int cwipiStep = readInt(runTime.controlDict().lookup("cwipiStep"));
int cwipiTimeStep = cwipiStep;

volPointInterpolation pInterp(mesh);

volVectorField momentumEquation(
    IOobject(
        "momentumEquation",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE),
    fvc::curl(U) ^ U);
// volScalarField continuityEquation(
//     IOobject(
//         "continuityEquation",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE),
//     ((thermo.rho() * Foam::pow(c, 2)) / thermo.Cp()) * (fvc::ddt(s) + (U & fvc::grad(s))));

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
