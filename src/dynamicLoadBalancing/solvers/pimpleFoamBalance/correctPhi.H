DebugVar("LLLL");

/*
refPtr<surfaceScalarField> tmeshPhi = mesh.setPhi();

if (tmeshPhi)
{
    DebugVar("MESHPHI EXISTS");
}
else
{
    DebugVar("MESHPHI DOES NOT EXIST");

    tmeshPhi.reset
    (
        std::make_unique<surfaceScalarField>
        (
        IOobject
        (
            "meshPhi",
            mesh.time().timeName(),
            mesh,
            IOobject::LAZY_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimVolume/dimTime, Foam::zero{})
        )
    );

    // tmeshPhi.ref() = Zero;

    DebugVar("JJJJ");

    if (tmeshPhi)
    {
        DebugVar("NOW MESHPHI EXISTS");
    }
    else
    {
        DebugVar("STILL MESHPI DOES NOT EXIST");
    }
}
*/
CorrectPhi
(
    U,
    phi,
    p,
    dimensionedScalar("rAUf", dimTime, 1),
    geometricZeroField(),
    pimple
);

DebugVar("MMMM");

#include "continuityErrs.H"

DebugVar("NNNN");
