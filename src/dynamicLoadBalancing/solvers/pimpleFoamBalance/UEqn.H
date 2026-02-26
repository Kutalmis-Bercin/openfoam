// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) + fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}

/*
        if (UPstream::myProcNo() == 2)
        {
            volatile double d = 0;
            for (int n = 0; n != 10000; ++n)
            {
                const scalar e = rng.sample01<scalar>();
                for (int m = 0; m != 100000; ++m)
                {
                    d += d * n * m * e;
                }
            }
        }
*/