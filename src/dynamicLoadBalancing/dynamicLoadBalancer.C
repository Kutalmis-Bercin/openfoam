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

#include "dynamicLoadBalancer.H"
#include "balanceModel.H"
#include "performanceMetric.H"
#include "profilingPstream.H"
#include "Function1.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicLoadBalancer, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::dynamicLoadBalancer::isBalancingDue() const
{

    Pout<< "# dynamicLoadBalancer::isBalancingDue(): Checking if balancing is due..."
        << endl;


    // Requires that  isImbalanced() has been called before

    const scalar t = mesh_.time().value();
    const label timeIndex = mesh_.time().timeIndex();
    const label freq = max(label(1), frequencyPtr_().value(t));

    bool skipBalance = false;
    if (skipBalancePtr_)
    {
        skipBalance = bool(skipBalancePtr_().value(t));
    }


    Pout<< "# Current time index: " << timeIndex << nl
        << "# time value: " << t << nl
        << "# Balancing frequency: " << freq << nl;


    bool isBalancingDue =
        (
            UPstream::nProcs() > 1
         && timeIndex > 1
         && curTimeIndex_ != timeIndex
         && (freq <= 1 || (timeIndex % freq == 0))
         && !skipBalance
        );

    if (isBalancingDue)
    {
        if (checkChronicImbalance_)
        {
            const scalar tol = recomputeImbalanceTolerancePtr_();
            const scalar backoffPeriod = backoffPeriodPtr_();

            // Already reduced across processors
            const scalar currImbl = performanceMetricPtr_->imbalanceValue();
            const scalar prevImbl = performanceMetricPtr_->prevImbalanceValue();

            if (mag(currImbl - prevImbl) < tol)  // this clashes with imbalanceThreshold
            {
                if ((t - lastBalanceTime_) < backoffPeriod)
                {
                    Info<< "    Chronic imbalance: skipping balancing" << endl;

                    isBalancingDue = false;
                }
            }
        }
    }

    return isBalancingDue;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicLoadBalancer::dynamicLoadBalancer
(
    fvMesh& mesh
)
:
    mesh_(mesh),
    dict_
    (
        IOdictionary::readContents
        (
            IOobject
            (
                "dynamicLoadBalanceDict",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ
            )
        )
    ),
    balanceModelPtr_(balanceModel::New(mesh, dict_)),
    performanceMetricPtr_(performanceMetric::New(mesh, dict_)),
    frequencyPtr_(nullptr),
    skipBalancePtr_(nullptr),
    recomputeImbalanceTolerancePtr_(nullptr),
    backoffPeriodPtr_(nullptr),
    checkChronicImbalance_(false),
    lastBalanceTime_(0.0),
    curTimeIndex_(-1)
{

    Pout<< "dynamicLoadBalancer::dynamicLoadBalancer() : Enabling "
        << "MPI profiling for dynamic load balancing." << endl;


    profilingPstream::enable();
    read(dict_);


    Pout<< "dynamicLoadBalancer::dynamicLoadBalancer() : Dynamic load balancer "
        << "constructed successfully." << endl;


    if (!UPstream::parRun())
    {
        WarningInFunction
            << "Dynamic load balancing requires running in parallel."
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicLoadBalancer::~dynamicLoadBalancer()
{
    profilingPstream::disable();
}  // Various models were forward declared


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicLoadBalancer::read(const dictionary& dict)
{

    Pout<< "dynamicLoadBalancer::read() : Reading dynamic load balancer "
        << " dictionary" << endl;

    if
    (
        !balanceModelPtr_->read(dict)
     || !performanceMetricPtr_->read(dict)
    )
    {
        return false;
    }

    skipBalancePtr_ =
        Function1<label>::NewIfPresent("skipBalance", dict, &mesh_);


    if
    (
        const auto& chronic = dict.subOrEmptyDict("chronicImbalanceSettings");
       !chronic.empty()
    )
    {
        recomputeImbalanceTolerancePtr_ = autoPtr<scalar>::New
        (
            chronic.getScalar("recomputeImbalanceTolerance")
        );
        backoffPeriodPtr_ = autoPtr<scalar>::New
        (
            chronic.getScalar("backoffPeriod")
        );

        if (recomputeImbalanceTolerancePtr_ && backoffPeriodPtr_)
        {
            checkChronicImbalance_ = true;
        }
    }



    Pout<< "dynamicLoadBalancer::read() : Setting balancing frequency"
        << endl;


    frequencyPtr_ = Function1<label>::New("frequency", dict, &mesh_);

    bool ok = true;
    ok = ok && (frequencyPtr_);

    return ok;
}


bool Foam::dynamicLoadBalancer::balance()
{

//
    Pout<< "dynamicLoadBalancer::balance() : Checking if balancing is due..." << endl;
//

    const bool imbalanced = performanceMetricPtr_->isImbalanced();
    const bool due = isBalancingDue();
    const bool anyDueImbalanced = returnReduce(due && imbalanced, orOp<bool>());

//
    Pout<< tab << "imbalanced: " << imbalanced << endl;
    Pout<< tab << "due: " << due << endl;
    Pout<< tab << "anyDueImbalanced: " << anyDueImbalanced << endl;
//

/*
    if (!anyDueImbalanced)
    {
        return false;
    }
*/
    curTimeIndex_ = mesh_.time().timeIndex();
    lastBalanceTime_ = mesh_.time().value();

//
    Pout<< tab << "# Dynamic load balancing is going to be performed" << endl;
//

    //mesh_.topoChanging(false);

    bool balanced =
        balanceModelPtr_->balance(performanceMetricPtr_->processorWeights());

    Pout<< tab << "# DLB: DONE" << endl;

    performanceMetricPtr_->resetTimeCounters();

    Pout<< tab<< "# dynamicLoadBalancer: resetTimeCounters was reset." << endl;

    return balanced;
}


// ************************************************************************* //
