if (cwipiTimeStep == cwipiStep)
{
    Info << "Before wait." << endl;
    cwipiWait();
    Info << "After wait." << endl;
    cwipiTimeStep = 0;

    // Gradually ramp up to source strength of 100%
    // if (cwipiRampIter < cwipiRampSteps)
    // {
    //     cwipiRampIter = cwipiRampIter + 1;
    //     cwipiRampStrength = static_cast<Foam::scalar>(cwipiRampIter) / static_cast<Foam::scalar>(cwipiRampSteps);
    // }
    // else
    // {
    //     cwipiRampIter = cwipiRampSteps;
    //     cwipiRampStrength = 1.0;
    // }

    // Info << "CWIPI source strength = " << cwipiRampStrength << endl;
}
cwipiTimeStep = cwipiTimeStep + 1;
