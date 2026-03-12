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

#include "performanceMetric.H"
#include "fvMesh.H"
#include "Function1.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(performanceMetric, 0);
    defineRunTimeSelectionTable(performanceMetric, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::performanceMetric::performanceMetric
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    imbalanceThresholdPtr_(nullptr),
    minUsage_(1e-3),
    smoothing_(-1),
    localUsage_(0),
    imbalanceValue_(0),
    prevImbalanceValue_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::performanceMetric::~performanceMetric()
{}  // fvMesh is forward declared


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::performanceMetric::read(const dictionary& dict)
{
#if DEBUGME
    Pout<< "performanceMetric::read() : Reading performance metric "
        << " dictionary" << endl;
#endif

    smoothing_ = dict.getOrDefault<scalar>("smoothing", -1);
    minUsage_ = dict.getOrDefault<scalar>("minUsage", minUsage_);
    minUsage_ = max(minUsage_, SMALL);


    imbalanceThresholdPtr_ =
        Function1<scalar>::New("imbalanceThreshold", dict, &mesh_);

    bool ok = true;
    ok = ok && imbalanceThresholdPtr_;

    return ok;
}


Foam::scalar Foam::performanceMetric::calcImbalance() const
{
    // Update local usage from the current mesh state
    localUsage_ = this->sampleLocalUsage();

    const label nProcs = UPstream::nProcs();

    // Sum usage across ranks
    const scalar sumUsage = returnReduce(localUsage_, sumOp<scalar>());

    if (nProcs <= 1 || sumUsage <= VSMALL)
    {
        return 0.0;
    }

    // Arithmetic average usage
    const scalar avg = sumUsage/scalar(nProcs);

    // Local relative deviation from the average
    const scalar localDev = mag(1.0 - localUsage_/(avg + VSMALL));

    // Global maximum of the relative deviation
    const scalar rawMax = returnReduce(localDev, maxOp<scalar>());

    // Theoretical maximum is (nProcs-1), scale to [0,1]
    const scalar denom = scalar(nProcs - 1);
    const scalar scaledMax =
        Foam::min
        (
            Foam::max
            (
                rawMax/(denom + VSMALL),
                scalar(0)
            ),
            scalar(1)
        );

    return scaledMax;
}


bool Foam::performanceMetric::isImbalanced() const
{
#if DEBUGME
    Pout<< "performanceMetric::isImbalanced() : Checking imbalance..." << endl;
#endif

    imbalanceValue_ = calcImbalance();

    const scalar t = mesh_.time().value();
    const scalar imbalanceThreshold = imbalanceThresholdPtr_().value(t);

#if DEBUGME
    Pout<< "Imbalance threshold: " << imbalanceThreshold << endl;
#endif

#if DEBUGME
    Pout<< "Calculated imbalance value: " << imbalanceValue_ << endl;
#endif

    bool isImbalanced = (imbalanceValue_ > imbalanceThreshold);

    // Already reduced imbalance values - MOVED?
    prevImbalanceValue_ = imbalanceValue_;

    return isImbalanced;
}


Foam::tmp<Foam::scalarField> Foam::performanceMetric::processorWeights() const
{
    const label nProcs = UPstream::nProcs();
    scalarField usage(nProcs, 0);

    for (label i = 0; i < nProcs; ++i)
    {
        if (UPstream::myProcNo() == i)
        {
            usage[i] = localUsage_;
        }
    }
    Pstream::listReduce(usage, maxOp<scalar>());


    // Inverse-usage capacity and normalise
    auto tweights = tmp<scalarField>::New(nProcs, scalar(0));
    scalarField& weights = tweights.ref();
    scalar sumCap = 0;
    forAll(weights, i)
    {
        const scalar cap = 1.0/max(minUsage_, usage[i]);
        weights[i] = cap;
        sumCap += cap;
    }


    if (sumCap <= SMALL)
    {
        // Fallback to uniform target if all zero
        weights = scalarField(nProcs, scalar(1.0/scalar(max(1, nProcs))));
    }
    else
    {
        // Larger weights for higher usage
        weights /= sumCap;
    }

    // Already reduced imbalance values - MOVED?
    // prevImbalanceValue_ = imbalanceValue_;


    auto tweightsTest = tmp<scalarField>::New(nProcs, scalar(0));
    scalarField& weightsTest = tweightsTest.ref();
    weightsTest[0] = 998959;  // problematic processor weights
    weightsTest[1] = 1041;

    return tweightsTest;
//    return tweights;
}


// ************************************************************************* //
