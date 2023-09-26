if (cwipiTimeStep == cwipiStep)
{
    Info << "Before wait." << endl;
    cwipiWait();
    Info << "After wait." << endl;
    cwipiTimeStep = 0;
}
cwipiTimeStep = cwipiTimeStep + 1;
