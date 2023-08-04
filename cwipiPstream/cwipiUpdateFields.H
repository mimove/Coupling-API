L = fvc::curl(U) ^ U;
s = (thermo.Cp() * log(thermo.T() / Foam::constant::standard::Tstd)) - ((Foam::constant::physicoChemical::R / thermo.W()) * log(thermo.p() / Foam::constant::standard::Pstd));
c = Foam::sqrt(thermo.gamma() * thermo.p() / thermo.rho());
