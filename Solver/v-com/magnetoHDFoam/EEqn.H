{
    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he) 
      + fvm::div(phi, he)
      - fvm::laplacian(turbulence->alphaEff(), he)
      - ETAMU*magSqr(fvc::curl(B))
    ==
        fvOptions(rho, he)
     );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();
}
