{
    volScalarField viscDiss
    (
        "viscDiss",
        turbulence->devRhoReff(U) && fvc::grad(U) // tau:grad(U)
    );

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, e) 
      + fvm::div(phi, e)
      - fvm::laplacian(turbulence->alphaEff(), e)
     ==
      - p * fvc::div(U)                 
      + viscDiss                          
      + ETAMU * magSqr(fvc::curl(B))
      + fvOptions(rho, e)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(e);

    thermo.correct();
}

