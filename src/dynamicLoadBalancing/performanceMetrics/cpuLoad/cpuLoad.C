/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "cpuLoad.H"
#include "Pstream.H"
#include "profilingPstream.H"
#include <algorithm>
#include <limits>
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace performanceMetrics
{
    defineTypeNameAndDebug(cpuLoad, 0);
    addToRunTimeSelectionTable
    (
        performanceMetric,
        cpuLoad,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::performanceMetrics::cpuLoad::cpuLoad
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    performanceMetric(mesh, dict),
    clock_(),
    cpuTime_()
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::performanceMetrics::cpuLoad::read(const dictionary& dict)
{
    if (!performanceMetric::read(dict))
    {
        return false;
    }

    return true;
}


Foam::scalar Foam::performanceMetrics::cpuLoad::sampleLocalUsage() const
{
    Pout<< "performanceMetrics::cpuLoad::sampleLocalUsage() : "
        << "Sampling local CPU load..." << endl;
    // Measure wall-clock time increment since the previous restart
    scalar dLocalWall = clock_.clockTimeIncrement();

    Pout<< "Local clock time increment: " << dLocalWall << nl;

    // A gate where all processors wait the slowest one
    scalar dWall = returnReduce(dLocalWall, maxOp<scalar>());

    Pout<< "Global wall time increment: " << dWall << nl;

    // Measure the simplistic MPI profiling times after all the processors
    // have reached this point thanks to the returnReduce above
    scalar dMpi = 0;
    const profilingPstream::timingList& mpiTimes = profilingPstream::times();

    // Pout<< "MPI time increments: " << dMpi << nl;

    for (unsigned i = 0; i < profilingPstream::nCategories; ++i)
    {
        dMpi += mpiTimes[i];
    }

    Pout<< "Total MPI time increment: " << dMpi << nl;

    // Measure CPU time increment since the previous restart; ideally, the
    // dCpu includes the MPI profiling times as well
    const scalar dCpuAndMpi  = cpuTime_.cpuTimeIncrement();

    // Pout<< "CPU time & MPI time increment: " << dCpuAndMpi << nl;

    // Estimate the CPU time increment excluding the MPI times
    scalar dCpu = max(dCpuAndMpi - dMpi, 0.0);

    Pout<< "CPU time increment: " << dCpu << nl;

    // Avoid problems with zero wall-clock time increments
    if (dWall <= std::numeric_limits<scalar>::epsilon())
    {
        dWall = 1.0;
        dCpu  = 0.0;
    }

    Pout<< "dWall adjusted: " << dWall << nl
        << "dCpu adjusted: " << dCpu << nl;

    // Warn the user if dCpu is larger than dWall
#ifdef FOAM_DEBUG
    if (dCpu > dWall)
    {
        WarningInFunction
            << "Estimated CPU time increment (" << dCpu << " s) is larger "
            << "than wall-clock time increment (" << dWall << " s). "
            << "This may lead to CPU load values larger than 1.0."
            << endl;
    }
#endif

    scalar usage = dCpu/dWall;

    Pout<< "Raw local usage (dCpu/dWall): " << usage << nl;

    // Smooth usage to avoid large jumps
    if (smoothing() > SMALL && smoothing() <= 1)
    {
        static scalar prevUsage = usage;
        usage = smoothing()*usage + (1 - smoothing())*prevUsage;
        prevUsage = usage;
    }

    Pout<< "Clock time increment: " << dLocalWall << nl
        << "Global wall time increment: " << dWall << nl
        << "CPU time & MPI time increment: " << dCpuAndMpi << nl
        << "MPI time increment: " << dMpi << nl
        << "CPU time increment: " << dCpu << nl
        << "Local usage: " << usage << nl;

    // Clamp to minimum
    usage = max(usage, scalar(minUsage()));
    return usage;
}


bool Foam::performanceMetrics::cpuLoad::resetTimeCounters()
{
    // The 'clock' is reset internally in clockTimeIncrement()

    cpuTime_.resetCpuTimeIncrement();

    profilingPstream::reset();

    return true;
}


// ************************************************************************* //
