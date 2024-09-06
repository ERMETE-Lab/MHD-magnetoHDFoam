fvVectorMatrix BEqn
(
    fvm::ddt(B)
    + B*fvc::div(U)
    + fvm::div(linearInterpolate(U) & mesh.Sf(), B)
    - fvc::div(phiB, U)
    - fvm::laplacian(DB, B)
);

BEqn.solve();

volScalarField rAB(1.0/BEqn.A());
surfaceScalarField rABf("rABf", fvc::interpolate(rAB));

phiB = fvc::flux(B);

while (bpiso.correctNonOrthogonal())
{
    fvScalarMatrix pBEqn
    (
        fvm::laplacian(rABf, pB) == fvc::div(phiB)
    );

    pBEqn.solve();

    if (bpiso.finalNonOrthogonalIter())
    {
        phiB -= pBEqn.flux();
    }
}
