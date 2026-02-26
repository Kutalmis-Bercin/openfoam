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

        // To test dynamic load balancing cpuLoad performance metric
/*
        if (UPstream::myProcNo() == 0)
        {
            volatile double d = 0;
            // for (int n = 0; n != 10000; ++n)
            for (int n = 0; n != 100; ++n)
            {
                const scalar e = rng.sample01<scalar>();
                for (int m = 0; m != 100000; ++m)
                {
                    d += d * n * m * e;
                }
            }
        }
*/

        /*
        Pout
            << "CPU time increment: " << runTime.cpuTimeIncrement() << nl
            << "Clock time increment: " << runTime.clockTimeIncrement() << nl;
        */

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
//std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();
                //DebugVar("pimpleFoam: Entering loadBalancer");

                //DebugVar("pimpleFoam: Exiting loadBalancer");
//std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
//std::cout << "Time difference (sec) = " <<  (std::chrono::duration_cast<std::chrono::microseconds>(end1 - begin1).count()) /1000000.0  <<std::endl;

                DebugVar("pimpleFoam: before mesh.changing");
                if (mesh.changing())
                {
                    DebugVar("pimpleFoam: mesh.changing");
                    MRF.update();

                    if (correctPhi)
                    {
                        DebugVar("pimpleFoam: correctPhi");
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);

                        DebugVar("pimpleFoam: makeRelative");
                    }

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }

                    mesh.topoChanging(false);
                }

                DebugVar("pimpleFoam: before loadBalancer");

                // mesh.clearMeshPhi();
                loadBalancer.balance();

            }

//std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
            #include "UEqn.H"
//std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
//std::cout << "Time difference (sec) = " <<  (std::chrono::duration_cast<std::chrono::microseconds>(end2 - begin2).count()) /1000000.0  <<std::endl;

            // --- Pressure corrector loop
//std::chrono::steady_clock::time_point begin3 = std::chrono::steady_clock::now();
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
//std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();
//std::cout << "Time difference (sec) = " <<  (std::chrono::duration_cast<std::chrono::microseconds>(end3 - begin3).count()) /1000000.0  <<std::endl;

            if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
