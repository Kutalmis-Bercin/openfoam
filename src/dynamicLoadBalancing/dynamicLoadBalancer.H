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
    Foam::dynamicLoadBalancer

Description
    Performs the dynamic load balancing of the mesh and fields according to the
    specified performance-metric, load and balance models.

Usage
    Minimal example by using \c constant/dynamicMeshDict:
    \verbatim
    dynamicLoadBalancerCoeffs
    {
        // Mandatory entries
        balanceModel        <word>;
        performanceMetric   <word>;
        frequency           <Function1<label>>;

        // Optional entries
        skipBalance         <Function1<label>>;
        checkChronicImbalance <Function1<label>>;

        // Inherited entries
        imbalanceThreshold  <Function1<scalar>>;
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                        | Type | Reqd | Deflt
      balanceModel | Name of specified balance model    | word | yes  | -
      performanceMetric | Name of specified performance metric | word | yes | -
      frequency | Frequency (in time step/iteration index) when      <!--
            --> the dynamic load balancing is performed | F1<label> | yes | -
      imbalanceThreshold | Threshold value that defines when load balancing <!--
                --> should be triggered              | F1<scalar> | yes | -
      skipBalance | Skip load balancing if true      | F1<label>  | no  | false
      checkChronicImbalance | Automatically skip load balancing if true
                -->                                  | F1<label>  | no  | false
    \endtable

    Options for the \c balanceModel entry:
    \verbatim
      distributed  | Distributes the mesh in a decentralised manner
    \endverbatim

    Options for the \c performanceMetric entry:
    \verbatim
      cellNumber   | Identifies imbalance based on cell-number distribution
      cpuLoad      | Identifies imbalance based on CPU load distribution
    \endverbatim

Note
  - The dynamic load balancer requires running in parallel.
  - Limited for the finite-volume meshes only.
  - The \c frequency variable ignore values less than 1 and treats them as 1.

SourceFiles
    dynamicLoadBalancer.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_dynamicLoadBalancer_H
#define Foam_dynamicLoadBalancer_H

#include "fvMesh.H"

#define DEBUGME 1

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class balanceModel;
class performanceMetric;
template <class Type> class Function1;

/*---------------------------------------------------------------------------*\
                    Class dynamicLoadBalancer Declaration
\*---------------------------------------------------------------------------*/

class dynamicLoadBalancer
{
    // Private Data

        //- Reference to the finite-volume mesh
        fvMesh& mesh_;

        //- Top-level dictionary for the user settings
        dictionary dict_;

        //- Load balance model that manages the load balance across processors
        //- by redistributing the mesh and fields on request
        autoPtr<balanceModel> balanceModelPtr_;

        //- Performance metric that defines the criteria for load balancing
        autoPtr<performanceMetric> performanceMetricPtr_;

        //- Frequency (time step/iteration index) at which the load balancing is
        //- performed
        autoPtr<Function1<label>> frequencyPtr_;

        //- Skip load balancing if true (0=false, non-zero=true)
        //  For chronic imbalance situations, repeated load balancing can be
        //  counterproductive; this option allows to skip balancing even though
        //  the imbalance criteria are met
        autoPtr<Function1<label>> skipBalancePtr_;

        //- Tolerance for recomputing imbalance for chronic imbalance
        autoPtr<scalar> recomputeImbalanceTolerancePtr_;

        //- Backoff period for load balancing attempts for chronic imbalance [s]
        autoPtr<scalar> backoffPeriodPtr_;

        //- Flag to indicate if chronic imbalance checking is enabled
        bool checkChronicImbalance_;

        //- The time of the last successful load balance attempt
        scalar lastBalanceTime_;

        //- The current time-step/iteration index
        label curTimeIndex_;


    // Private Member Functions

        //- Return true if the dynamic load balancing should be triggered
        bool isBalancingDue() const;


public:

    //- Declare type-name (with debug switch)
    ClassName("dynamicLoadBalancer");


    // Constructors

        //- Construct from the mesh
        dynamicLoadBalancer(fvMesh& mesh);


    // Generated Methods

        //- No copy construct
        dynamicLoadBalancer(const dynamicLoadBalancer&) = delete;

        //- No copy assignment
        void operator=(const dynamicLoadBalancer&) = delete;


    //- Destructor
    ~dynamicLoadBalancer();


    // Member Functions

    // Access

        //- Return reference to the top-level user dictionary
        const dictionary& dict() const noexcept { return dict_; }


    // Evaluation

        //- Perform the dynamic load balancing and return the status
        bool balance();


    // I/O

        //- Read the top-level dictionary
        bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
