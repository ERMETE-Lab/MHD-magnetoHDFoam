dimensionedScalar compressibility = fvc::domainIntegrate(psi);
bool compressible = (compressibility.value() > SMALL);

if (!pimple.SIMPLErho())
{
    rho = thermo.rho();
}

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
const volScalarField psip0(psi*p);

volScalarField rAU(1.0/UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho*rAU));
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p_rgh));

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phig(-rhorAUf*ghf*fvc::snGrad(rho)*mesh.magSf());


surfaceScalarField phiHbyA
(
    "phiHbyA",
    fvc::interpolate(rho)*fvc::flux(HbyA)
  + MRF.zeroFilter(rhorAUf*fvc::ddtCorr(rho, U, phi, rhoUf))
  + phig
);

fvc::makeRelative(phiHbyA, rho, U);
MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

// Update the pressure BCs to ensure flux consistency
constrainPressure(p_rgh, rho, U, phiHbyA, rhorAUf, MRF);

fvScalarMatrix p_rghDDtEqn
(
    fvc::ddt(rho) + psi*correction(fvm::ddt(p_rgh))
  + fvc::div(phiHbyA)
 ==
    fvOptions(psi, p_rgh, rho.name())
);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix p_rghEqn(p_rghDDtEqn - fvm::laplacian(rhorAUf, p_rgh));

    p_rghEqn.solve(p_rgh.select(pimple.finalInnerIter()));

    if (pimple.finalNonOrthogonalIter())
    {
        phi = phiHbyA + p_rghEqn.flux();
            
        // Explicitly relax pressure for momentum corrector
        p_rgh.relax();

        // Correct the momentum source with the pressure gradient flux
        // calculated from the relaxed pressure
        U = HbyA + rAU*fvc::reconstruct((phig + p_rghEqn.flux())/rhorAUf);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
        K = 0.5*magSqr(U);
    }
}

p = p_rgh + rho*gh;

pressureControl.limit(p);

if (p_rgh.needReference())
{
    if (!compressible)
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
    }
    else
    {
        p += (initialMass - fvc::domainIntegrate(psi*p))
            /compressibility;
        thermo.correctRho(psi*p - psip0, rhoMin, rhoMax);
        rho = thermo.rho();
        p_rgh = p - rho*gh;
        p_rgh.correctBoundaryConditions();
    }
}
else
{
    thermo.correctRho(psi*p - psip0, rhoMin, rhoMax);
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

rho = thermo.rho();

// Correct rhoUf if the mesh is moving
fvc::correctRhoUf(rhoUf, rho, U, phi);

if (thermo.dpdt())
{
    dpdt = fvc::ddt(p);

    if (mesh.moving())
    {
        dpdt -= fvc::div(fvc::meshPhi(rho, U), p);
    }
}
