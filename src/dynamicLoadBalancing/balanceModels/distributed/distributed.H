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
    Foam::balanceModels::distributed

Description
    The \c distributed balance model balances the mesh and field distribution
    by computing a new decomposition distribution from the supplied cell weights
    and triggering the master-processor distributor to send and receive the
    mesh and field parts.

Usage
    Minimal example by using \c constant/dynamicMeshDict:
    \verbatim
    dynamicLoadBalancerCoeffs
    {
        // Inherited entries
        ...

        // Mandatory entries
        balanceModel  distributed;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                            | Type | Reqd | Deflt
      type      | Model name: distributed                | word | yes  | -
    \endtable

SourceFiles
    distributed.C

\*---------------------------------------------------------------------------*/

#ifndef balanceModels_distributed_H
#define balanceModels_distributed_H

#include "balanceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace balanceModels
{

/*---------------------------------------------------------------------------*\
                            Class distributed Declaration
\*---------------------------------------------------------------------------*/

class distributed
:
    public balanceModel
{
    // Private Data

        //- Local dictionary for run-time distribution parameters
        autoPtr<IOdictionary> distributionDictPtr_;


    // Private Member Functions

        //- Return a labelList which is a transformation of a scalarList
        //- approximating the original scalar ratios as closely as possible.
        //  The 'processorWeights' entry only accepts label weights; therefore,
        //  the scalar CPU-load weights need to be transformed into labels.
        labelList transformToIntegers
        (
            const scalarField& w,
            const label targetSum = 1000000  // larger it is, smaller rel error
        );


public:

    //- Runtime type information
    TypeName("distributed");


    // Constructors

        //- Construct from components
        distributed
        (
            fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~distributed() = default;


    // Evaluation

        //- Balance the load by redistributing the mesh and fields
        virtual bool balance(const scalarField& weights);


    // I-O

        //- Read the top-level dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace balanceModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
