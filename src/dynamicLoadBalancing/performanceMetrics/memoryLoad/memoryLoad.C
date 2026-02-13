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

#include "memoryLoad.H"
#include "Pstream.H"
#include "profilingPstream.H"
#include "addToRunTimeSelectionTable.H"

#include <fstream>
#include <sstream>
#include <string>
#include <limits>
#include <cmath>

#if defined(__linux__)
    #include <unistd.h>
#elif defined(__APPLE__)
    #include <mach/mach.h>
#elif defined(_WIN32)
    #ifndef NOMINMAX
    #define NOMINMAX
    #endif
    #include <windows.h>
    #include <psapi.h>
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace performanceMetrics
{
    defineTypeNameAndDebug(memoryLoad, 0);
    addToRunTimeSelectionTable
    (
        performanceMetric,
        memoryLoad,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::performanceMetrics::memoryLoad::memoryLoad
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    performanceMetric(mesh, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::performanceMetrics::memoryLoad::read(const dictionary& dict)
{
    if (!performanceMetric::read(dict))
    {
        return false;
    }

    return true;
}


Foam::scalar Foam::performanceMetrics::memoryLoad::sampleLocalUsage() const
{
    // Return the current resident set size in megabytes (MB).
    scalar usageMB = 0.0;

#if defined(__linux__)

    // Preferred: /proc/self/statm (resident pages * page size)
    {
        std::ifstream statm("/proc/self/statm");
        if (statm.good())
        {
            unsigned long sizePages = 0, residentPages = 0;
            statm >> sizePages >> residentPages;

            const long pageSize = ::sysconf(_SC_PAGESIZE); // bytes
            if (pageSize > 0 && residentPages > 0)
            {
                usageMB = scalar(residentPages) * scalar(pageSize)
                        / scalar(1024.0 * 1024.0);
            }
        }
    }

    // Fallback: parse VmRSS from /proc/self/status (kB)
    if (usageMB <= SMALL)
    {
        std::ifstream status("/proc/self/status");
        std::string line;
        while (std::getline(status, line))
        {
            if (line.rfind("VmRSS:", 0) == 0)
            {
                std::istringstream iss(line.substr(6));
                double kb = 0.0;
                iss >> kb; // value in kB
                usageMB = scalar(kb / 1024.0);
                break;
            }
        }
    }

#elif defined(__APPLE__)

    // Use Mach task_info() to query resident size
    {
        mach_task_basic_info info;
        mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;

        if (task_info
            (
                mach_task_self(),
                MACH_TASK_BASIC_INFO,
                reinterpret_cast<task_info_t>(&info),
                &count
            ) == KERN_SUCCESS)
        {
            usageMB = scalar(info.resident_size) / scalar(1024.0 * 1024.0);
        }
    }

#elif defined(_WIN32)

    // Use PSAPI to query the working set size of the current process
    {
        PROCESS_MEMORY_COUNTERS pmc;
        if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
        {
            usageMB = scalar(pmc.WorkingSetSize) / scalar(1024.0 * 1024.0);
        }
    }

#else
    WarningInFunction
        << "memoryLoad: Unsupported platform for memory measurement. "
        << "Returning a constant usage of 1." << nl;
    usageMB = 1.0;
#endif

    // Guard against zeros and NaNs
    if (!std::isfinite(usageMB) || usageMB <= SMALL)
    {
        usageMB = 1.0; // minimal positive usage
    }

    return usageMB;
}


bool Foam::performanceMetrics::memoryLoad::resetTimeCounters()
{
    // ? profilingPstream::reset();

    return true;
}


// ************************************************************************* //
