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
    Foam::performanceMetrics::combinedLoad

Description
    A composite performance-metric that combines several existing
    performance-metric submodels (eg, \c cellNumber, \c cpuLoad) using
    user-specified weights. The composed metric evaluates the local usage
    as a weighted sum of the submodels' local usages

        usage = sum_i ( w_i * usage_i )

    and then leverages the base-class logic to compute per-processor target
    weights and imbalance. Submodels remain independent and are not merged.

Usage
    In the dynamic load-balancer dictionary:

    \verbatim
    performanceMetric  weightedCombination;

    weightedCombinationCoeffs
    {
        // List of submodels to combine
        models
        (
            {
                performanceMetric  cellNumber;
                weight  0.50;     // optional, defaults to equal weights
                // cellNumberCoeffs { ... } // if needed by submodel
            }
            {
                performanceMetric  cpuLoad;
                weight  0.50;     // optional
                // cpuLoadCoeffs { ... }
            }
        );
    }
    \endverbatim

    The generic entries supported by the base class (eg, \c minUsage,
    \c smoothing, \c imbalanceThreshold ...) can be placed at the same level
    as \c performanceMetric and \c weightedCombinationCoeffs.

SourceFiles
    combinedLoad.C

\*---------------------------------------------------------------------------*/

#ifndef performanceMetrics_combinedLoad_H
#define performanceMetrics_combinedLoad_H

#include "performanceMetric.H"
#include "PtrList.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace performanceMetrics
{

/*---------------------------------------------------------------------------*\
                            Class combinedLoad Declaration
\*---------------------------------------------------------------------------*/

class combinedLoad
:
    public performanceMetric
{
    // Private Data

        //- Submodels
        PtrList<performanceMetric> models_;

        //- Normalized, non-negative weights (one per submodel)
        List<scalar> weights_;


    // Private Member Functions

        //- Clear submodels/weights
        void clear();

        //- Construct submodels/weights from dictionary
        void constructFromDict(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("combinedLoad");


    // Constructors

        //- Construct from mesh and dictionary
        combinedLoad
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~combinedLoad() = default;


    // Member Functions

    // Evaluation

        //- Sample local usage (updates previous samples if smoothing>0)
        scalar sampleLocalUsage() const;

        //- Return true after resetting counters (no-op for this metric)
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
