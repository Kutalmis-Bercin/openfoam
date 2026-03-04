/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pimpleFoam.C

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.

    \heading Solver details
    The solver uses the PIMPLE (merged PISO-SIMPLE) algorithm to solve the
    continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}} + \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

Note
   The motion frequency of this solver can be influenced by the presence
   of "updateControl" and "updateInterval" in the dynamicMeshDict.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "dynamicLoadBalancer.H"
#include "Random.H"

#include <chrono>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );

    Random rng(1234567);

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "CourantNo.H"
            #include "setDeltaT.H"
        }

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        if (correctPhi)
        {
            Info<< "correctPhi is ON" << endl;
        }
        else
        {
            Info<< "correctPhi is OFF" << endl;
        }

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                // Do any mesh changes
                DebugVar("pimpleFoam: before any mesh changes");
                mesh.controlledUpdate();

                DebugVar("pimpleFoam: before mesh.changing");
                if (mesh.changing())
                {
                    DebugVar("AAAA");

                    mesh.clearMeshPhi();

                    DebugVar("BBBB");

                    loadBalancer.balance();

                    DebugVar("CCCC");
                    DebugVar("pimpleFoam: mesh.changing");


                    DebugVar("FFFF");

                    MRF.update();

                    DebugVar("HHHH");

                    if (correctPhi)
                    {
                        DebugVar("pimpleFoam: correctPhi");
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        DebugVar("DDDD");

                        #include "correctPhi.H"

                        DebugVar("EEEE");

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        DebugVar("pimpleFoam: makeRelative");
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }

                DebugVar("pimpleFoam: before loadBalancer");
            }


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

                refPtr<surfaceScalarField> tmeshPhi2 = mesh.setPhi();
                if (tmeshPhi2)
                {
                    DebugVar("NOW MESHPHI EXISTS");
                }
                else
                {
                    DebugVar("STILL MESHPI DOES NOT EXIST");
                }
            }

            DebugVar("1111");
            #include "UEqn.H"
            DebugVar("2222");

            DebugVar("3333");
            while (pimple.correct())
            {
                DebugVar("4444");
                #include "pEqn.H"
                DebugVar("5555");
            }

            DebugVar("6666");
            if (pimple.turbCorr())
            {
                DebugVar("7777");
                laminarTransport.correct();
                DebugVar("8888");
                turbulence->correct();
                DebugVar("9999");
            }
        }

        DebugVar("1010");
        runTime.write();
        DebugVar("2020");

        runTime.printExecutionTime(Info);
        DebugVar("3030");
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
