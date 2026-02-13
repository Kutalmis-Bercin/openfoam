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

Class
    Foam::performanceMetrics::memoryLoad

Description

Usage

Note

SourceFiles
    memoryLoad.C

\*---------------------------------------------------------------------------*/

#ifndef performanceMetrics_memoryLoad_H
#define performanceMetrics_memoryLoad_H

#include "performanceMetric.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace performanceMetrics
{

/*---------------------------------------------------------------------------*\
                            Class memoryLoad Declaration
\*---------------------------------------------------------------------------*/

class memoryLoad
:
    public performanceMetric
{
public:

    //- Runtime type information
    TypeName("memoryLoad");


    // Constructors

        //- Construct from mesh and dictionary
        memoryLoad
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~memoryLoad() = default;


    // Member Functions

    // Evaluation

        //- Sample local usage (updates previous samples)
        scalar sampleLocalUsage() const;

        //- Return true after resetting time counters
        //  NOTE: This may need to be refactored if more complex actions are
        //  needed later
        virtual bool resetTimeCounters();


    // I/O

        //- Read the top-level dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace performanceMetrics
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
