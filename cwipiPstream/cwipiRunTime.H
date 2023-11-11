coupling.updateFields(U, thermo);
if (coupling.sendNow() == true)
{
    coupling.updateSources(pInterp, thermo);
    F_0_p = pInterp.interpolate(coupling.F_p());
    F_0_u = pInterp.interpolate(coupling.F_u());
    coupling.send(F_0_p, F_0_u, mesh);
}
