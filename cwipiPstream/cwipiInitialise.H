coupling.createCoupling(mesh);
volPointInterpolation pInterp(mesh);
pointScalarField F_0_p = pInterp.interpolate(coupling.F_p());
pointVectorField F_0_u = pInterp.interpolate(coupling.F_u());
