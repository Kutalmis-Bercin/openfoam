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

#include "combinedLoad.H"
#include "Pstream.H"
#include "profilingPstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace performanceMetrics
{
    defineTypeNameAndDebug(combinedLoad, 0);
    addToRunTimeSelectionTable
    (
        performanceMetric,
        combinedLoad,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::performanceMetrics::combinedLoad::clear()
{
    models_.clear();
    weights_.clear();
}


void Foam::performanceMetrics::combinedLoad::constructFromDict
(
    const dictionary& dict
)
{
    // Parse own coeffs sub-dictionary
    const word coeffsName("combinedLoadCoeffs");
    const dictionary coeffs(dict.subOrEmptyDict(coeffsName));

    // The 'models' entry: a list of dictionaries
    if (!coeffs.found("models"))
    {
        FatalIOErrorInFunction(dict)
            << "Missing required entry 'models' in " << coeffsName << nl
            << exit(FatalIOError);
    }

    List<dictionary> modelDicts(coeffs.lookup("models"));

    const label n = modelDicts.size();

    models_.setSize(n);
    weights_.setSize(n);

    scalar sumW = 0;

    forAll(modelDicts, i)
    {
        dictionary& sub = modelDicts[i];

        // Each sub-dictionary must contain its own performanceMetric entry
        if (!sub.found("performanceMetric"))
        {
            FatalIOErrorInFunction(dict)
                << "Entry #" << i << " in '" << coeffsName
                << "' is missing 'performanceMetric' keyword\n"
                << "Sub-dict was: " << sub
                << exit(FatalIOError);
        }

        // Weight (optional), default to 1/n for now, normalize later
        scalar wi = 1;
        sub.readIfPresent("weight", wi);
        wi = max(wi, scalar(0)); // ensure non-negative
        weights_[i] = wi;
        sumW += wi;

        // Construct the submodel using the usual factory
        autoPtr<performanceMetric> pm = performanceMetric::New(mesh(), sub);
        models_.set(i, pm);
    }

    // Normalize weights. If all zero -> use equal weights
    if (sumW <= SMALL)
    {
        const scalar eq = n ? scalar(1.0)/scalar(n) : 0;
        forAll(weights_, i) { weights_[i] = eq; }
    }
    else
    {
        forAll(weights_, i) { weights_[i] /= sumW; }
    }
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::performanceMetrics::combinedLoad::combinedLoad
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    performanceMetric(mesh, dict),
    models_(),
    weights_()
{
    read(dict);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::performanceMetrics::combinedLoad::read(const dictionary& dict)
{
    if (!performanceMetric::read(dict))
    {
        return false;
    }

    clear();
    constructFromDict(dict);

    return true;
}


Foam::scalar Foam::performanceMetrics::combinedLoad::sampleLocalUsage() const
{
    if (models_.empty())
    {
        return scalar(1); // Avoid division by zero upstream
    }

    // Weighted sum of submodels' usages
    scalar usage = 0;
    forAll(models_, i)
    {
        const scalar ui = models_[i].sampleLocalUsage();
        usage += weights_[i]*ui;
    }

    // Optional smoothing from base (to avoid large jumps)
    if (smoothing() > SMALL && smoothing() <= 1)
    {
        static scalar prevUsage = usage;
        usage = smoothing()*usage + (1 - smoothing())*prevUsage;
        prevUsage = usage;
    }

    // Clamp with the base minimum usage
    usage = max(usage, scalar(minUsage()));

    return usage;
}


bool Foam::performanceMetrics::combinedLoad::resetTimeCounters()
{
    bool ok = true;
    forAll(models_, i)
    {
        ok = models_[i].resetTimeCounters() && ok;
    }
    return ok;
}


// ************************************************************************* //
