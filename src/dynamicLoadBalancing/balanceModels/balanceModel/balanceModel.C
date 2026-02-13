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

#include "balanceModel.H"
#include "decompositionMethod.H"
#include "decompositionModel.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(balanceModel, 0);
    defineRunTimeSelectionTable(balanceModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::balanceModel::balanceModel
(
    fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh)/*,
    decomposerPtr_
    (
        decompositionModel::New
        (
            mesh,
            IOdictionary
            (
                IOobject
                (
                    "decomposeParDict",
                    mesh_.time().system(),
                    mesh_.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                )
            )
        ).decomposer()
    )*/
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::balanceModel::~balanceModel()
{}  // fvMesh is forward declared


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::balanceModel::read(const dictionary& dict)
{
    return true;
}


// ************************************************************************* //
