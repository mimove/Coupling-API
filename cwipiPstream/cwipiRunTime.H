if (cwipiTimeStep == cwipiStep)
{
    LPrime = L - LMean;
    TPrime = T - TMean;
    sPrime = s - sMean;
    ThetaPrime = (TPrime * fvc::grad(sMean)) - (sPrime * fvc::grad(TMean));

    // Add Lamb vector if true
    if (lambVectorSwitch == true && entropyGradientSwitch == false)
    {
        momentumEquation = -LPrime * sourceDamping;
    }

    // Add entropy if true
    if (entropyGradientSwitch == true && lambVectorSwitch == false)
    {
        momentumEquation = ThetaPrime * sourceDamping;
    }

    // Add Lamb vector and entropy gradient if true
    if (lambVectorSwitch == true && entropyGradientSwitch == true)
    {
        momentumEquation = (ThetaPrime - LPrime) * sourceDamping;
    }

    acousticMomentumEquation = pInterp.interpolate(momentumEquation);

    if (entropyDerivativeSwitch == true)
    {
        continuityEquation = ((rhoMean / thermo.Cp()) * (fvc::ddt(sPrime) + (UMean & fvc::grad(sPrime)))) * sourceDamping;
        acousticContinuityEquation = pInterp.interpolate(continuityEquation);
    }

    cwipiReshapeSourceArrays(
        mesh,
        sourceArray,
        acousticContinuityEquation,
        acousticMomentumEquation,
        cwipiDim);

    cwipiSend(
        sourceArray,
        cwipiArgumentList,
        cwipiDim);
}
