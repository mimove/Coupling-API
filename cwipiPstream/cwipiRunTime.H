if (cwipiTimeStep == cwipiStep)
{
    LPrime = L - LMean;
    TPrime = T - TMean;
    sPrime = s - sMean;
    ThetaPrime = (TPrime * fvc::grad(sMean)) - (sPrime * fvc::grad(TMean));

    // Info << "CWIPI source strength = " << cwipiRampStrength << endl;

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

    // acousticMomentumEquation = pInterp.interpolate(cwipiRampStrength * momentumEquation);
    acousticMomentumEquation = pInterp.interpolate(momentumEquation);

    // if (entropyDerivativeSwitch == true)
    // {
    //     sPrime = s - sMean;
    //     continuityEquation = (((rhoMean * Foam::pow(cMean, 2)) / thermo.Cp()) * (fvc::ddt(sPrime) + (UMean & fvc::grad(sPrime)))) * sourceDamping;
    //     acousticContinuityEquation = pInterp.interpolate(continuityEquation);
    // }

    cwipiReshapeSourceArrays(mesh, sourceArray, acousticMomentumEquation, cwipiDim);

    cwipiSend(sourceArray, cwipiArgumentList, cwipiDim);
}
