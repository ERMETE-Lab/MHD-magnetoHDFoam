const volScalarField& psi = thermo.psi();
volScalarField& e = thermo.he(); // Energy balance in terms of internal energy e or he, it depends what you put in thermophysicalProperties in energy (thermotype)
