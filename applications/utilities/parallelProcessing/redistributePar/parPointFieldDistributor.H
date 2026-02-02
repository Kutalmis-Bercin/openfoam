/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2026 OpenCFD Ltd.
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
    Foam::parPointFieldDistributor

Description
    Distributor/redistributor for point fields,
    uses a two (or three) stage construction.

    The inconvenient multi-stage construction is needed since the
    pointMesh is directly associated with a polyMesh, which will probably
    have changed while creating the target mesh. This means that it is
    necessary to save the size of the source mesh and all of its
    patch meshPoints prior to making any changes (eg, creating the target
    mesh).

    -# Create with specified source mesh
    -# Save the meshPoints (per boundary) for the source mesh
    -# Attach a target mesh and mesh distribution
    -# Map the point fields
    .

    Runs in parallel. Redistributes from srcMesh to tgtMesh.

SourceFiles
    parPointFieldDistributor.cxx
    parPointFieldDistributor.txx

\*---------------------------------------------------------------------------*/

#ifndef Foam_parPointFieldDistributor_H
#define Foam_parPointFieldDistributor_H

#include "pointMesh.H"
#include "pointFieldsFwd.H"
#include "PtrList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapDistributePolyMesh;
class mapDistributeBase;
class IOobjectList;

/*---------------------------------------------------------------------------*\
                  Class parPointFieldDistributor Declaration
\*---------------------------------------------------------------------------*/

class parPointFieldDistributor
{
    // Private Data

        //- The source mesh reference
        const pointMesh& srcMesh_;

        //- Number of points in the old (source) mesh
        const label nOldPoints_;

        //- The pointPatch mesh points
        PtrList<labelList> patchMeshPoints_;

        //- The target (destination) mesh reference
        refPtr<pointMesh> tgtMeshRef_;

        //- Distribution map reference
        refPtr<mapDistributePolyMesh> distMapRef_;

        //- Point patch mappers
        PtrList<mapDistributeBase> patchPointMaps_;

        //- Storage for dummy handler (when using bool control)
        refPtr<fileOperation> dummyHandler_;

        //- Write control via a file handler
        refPtr<fileOperation>& writeHandler_;

        //- Write control as a bool
        Switch isWriteProc_;


public:

    //- Output verbosity when writing
    static int verbose_;


    // Generated Methods

        //- No copy construct
        parPointFieldDistributor(const parPointFieldDistributor&) = delete;

        //- No copy assignment
        void operator=(const parPointFieldDistributor&) = delete;


    // Constructors

        //- Basic construction
        explicit parPointFieldDistributor
        (
            //! The source pointMesh
            const pointMesh& srcMesh,
            //! Call saveMeshPoints() immediately. Normally false
            const bool savePoints,
            //! Tagged for output writing (on this proc)
            const bool isWriteProc
        );

        //- Basic construction
        explicit parPointFieldDistributor
        (
            //! The source pointMesh
            const pointMesh& srcMesh,
            //! Call saveMeshPoints() immediately. Normally false
            const bool savePoints,
            //! Valid for output writing (on this proc)
            refPtr<fileOperation>& writeHandler
        );

        //- Full construction of source/target
        explicit parPointFieldDistributor
        (
            //! The source pointMesh
            const pointMesh& srcMesh,
            //! The target pointMesh
            const pointMesh& tgtMesh,
            //! The distribution map
            const mapDistributePolyMesh& distMap,
            //! Call saveMeshPoints() immediately. Normally false
            const bool savePoints,
            //! Tagged for output writing (on this proc)
            const bool isWriteProc
        );

        //- Full construction of source/target
        explicit parPointFieldDistributor
        (
            //! The source pointMesh
            const pointMesh& srcMesh,
            //! The target pointMesh
            const pointMesh& tgtMesh,
            //! The distribution map
            const mapDistributePolyMesh& distMap,
            //! Call saveMeshPoints() immediately. Normally false
            const bool savePoints,
            //! Valid for output writing (on this proc)
            refPtr<fileOperation>& writeHandler
        );


    // Member Functions

        //- The source mesh
        const pointMesh& sourceMesh() const noexcept { return srcMesh_; }

        //- True if meshPoints (per boundary) for the source mesh
        //- have been saved
        bool hasMeshPoints() const noexcept;

        //- True if patch maps (per boundary) exist
        bool hasPatchPointMaps() const noexcept;

        //- True if a target mesh/distribution map has been attached
        bool hasTarget() const noexcept;

        //- Clear out meshPoints (per boundary) for the source mesh
        void clearMeshPoints();

        //- Clear out patch maps (per boundary)
        void clearPatchPointMaps();

        //- Create/recreate meshPoints (per boundary) for the source mesh
        void saveMeshPoints();

        //- Construct per-patch addressing
        void createPatchPointMaps();

        //- Clear target mesh / distribution map
        void resetTarget();

        //- Reset target mesh / distribution map
        void resetTarget
        (
            const pointMesh& tgtMesh,
            const mapDistributePolyMesh& distMap
        );


    // Field Mapping

        //- Read, distribute and write all/selected point field types
        //- (scalar, vector, ... types)
        label distributeAllFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        ) const;


        //- Distribute point field
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh>>
        distributeField
        (
            const GeometricField<Type, pointPatchField, pointMesh>& fld
        ) const;

        //- Read and distribute point field
        template<class Type>
        tmp<GeometricField<Type, pointPatchField, pointMesh>>
        distributePointField
        (
            const IOobject& fieldObject
        ) const;

        //- Read, distribute and write all/selected point fields
        template<class Type>
        label distributePointFields
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        ) const;

        //- Distributed each (unregistered!) point field
        //- and store the result on its objectRegistry
        template<class Type>
        void distributeAndStore
        (
            const PtrList<GeometricField<Type, pointPatchField, pointMesh>>&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "parPointFieldDistributor.txx"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
