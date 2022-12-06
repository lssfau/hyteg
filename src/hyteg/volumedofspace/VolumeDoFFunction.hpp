/*
* Copyright (c) 2017-2022 Nils Kohl.
*
* This file is part of HyTeG
* (see https://i10git.cs.fau.de/hyteg/hyteg).
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include "hyteg/ReferenceCounter.hpp"
#include "hyteg/celldofspace/CellDoFIndexing.hpp"
#include "hyteg/facedofspace_old/FaceDoFIndexing.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/volumedofspace/VolumeDoFIndexing.hpp"
#include "hyteg/volumedofspace/VolumeDoFPackInfo.hpp"

namespace hyteg {
namespace volumedofspace {

/// \brief Parallel data structure for degrees of freedom that live in volumes.
///
/// The VolumeDoFFunction has the following properties:
///  - a specified number of scalars is allocated per micro-volume
///  - in 2D DoFs are allocated on macro-faces only (micro-volumes are micro-faces)
///  - in 3D DoFs are allocated on macro-cells only (micro-volumes are micro-cells)
///  - ghost-layers are allocated so that each inner micro-volume of a macro-primitives has neighboring DoFs interfaced over
///    all facets (dim - 1, that is edges in 2D (3 facets) and faces in 3D (4 facets))
///  - the number of scalars allocated on each macro-volume can be specified individually to allow for the implementation of e.g.
///    adaptive p-refinement, the scalars on the ghost-layers are adapted accordingly
///  - there is no "position" assigned to the DoFs - meaning to the scalars has to be given by external data structures
template < typename ValueType >
class VolumeDoFFunction : public Function< VolumeDoFFunction< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = VolumeDoFFunction< VType >;

   /// \brief Allocates a VolumeDoFFunction.
   ///
   /// \param name                   name of this function
   /// \param storage                PrimitiveStorage object
   /// \param minLevel               smallest h-refinement level to allocate
   /// \param maxLevel               largest h-refinement level to allocate
   /// \param numScalarsPerPrimitive map of (local) PrimitiveIDs to the number of DoFs per micro-element that shall be allocated
   /// \param memoryLayout           Structure-of-Arrays or Array-of-Structures
   VolumeDoFFunction( const std::string&                         name,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      uint_t                                     minLevel,
                      uint_t                                     maxLevel,
                      const std::map< PrimitiveID, uint_t >&     numScalarsPerPrimitive,
                      indexing::VolumeDoFMemoryLayout            memoryLayout );

   /// \brief Allocates a VolumeDoFFunction.
   ///
   /// \param name         name of this function
   /// \param storage      PrimitiveStorage object
   /// \param minLevel     smallest h-refinement level to allocate
   /// \param maxLevel     largest h-refinement level to allocate
   /// \param numScalars   number of scalars allocated per element (same for all macro-volumes)
   /// \param memoryLayout Structure-of-Arrays or Array-of-Structures
   VolumeDoFFunction( const std::string&                         name,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      uint_t                                     minLevel,
                      uint_t                                     maxLevel,
                      uint_t                                     numScalars,
                      indexing::VolumeDoFMemoryLayout            memoryLayout );

   /// Destructor
   ~VolumeDoFFunction();

   /// Copy constructor
   VolumeDoFFunction( const VolumeDoFFunction< ValueType >& other );

   /// Copy assignment
   VolumeDoFFunction& operator=( const VolumeDoFFunction< ValueType >& other );

   /// \brief Updates ghost-layers.
   void communicate( uint_t level );

   /// \brief Assigns a linear combination of multiple VolumeDoFFunctions to this.
   void assign( const std::vector< ValueType >&                                                      scalars,
                const std::vector< std::reference_wrapper< const VolumeDoFFunction< ValueType > > >& functions,
                uint_t                                                                               level );

   /// \brief Evaluates the dot product on all local DoFs. No communication is involved and the results may be different on each
   /// process.
   ValueType dotLocal( const VolumeDoFFunction< ValueType >& rhs, uint_t level ) const;

   /// \brief Evaluates the (global) dot product. Involves communication and has to be called collectively.
   ValueType dotGlobal( const VolumeDoFFunction< ValueType >& rhs, uint_t level ) const;

   /// \brief Returns a pointer to the array that stores all degrees of freedom.
   ///
   /// You better know what you are doing when calling this ...
   ///
   /// \param primitiveID macro-primitive (volume primitive) to get the data of
   /// \return pointer to the array that stores the data
   ValueType* dofMemory( PrimitiveID primitiveID, uint_t level ) const
   {
      if ( this->storage_->hasGlobalCells() )
      {
         WALBERLA_CHECK( this->storage_->cellExistsLocally( primitiveID ),
                         "Cannot read/write DoF since macro-cell does not exists (locally)." );
         FunctionMemory< ValueType >* fmem = this->storage_->getCell( primitiveID )->getData( cellInnerDataID_ );
         WALBERLA_CHECK( fmem->hasLevel( level ), "Memory was not allocated for level " << level << "." );
         auto data = fmem->getPointer( level );
         return data;
      }
      else
      {
         WALBERLA_CHECK( this->storage_->faceExistsLocally( primitiveID ),
                         "Cannot read/write DoF since macro-face does not exists (locally)." );
         FunctionMemory< ValueType >* fmem = this->storage_->getFace( primitiveID )->getData( faceInnerDataID_ );
         WALBERLA_CHECK( fmem->hasLevel( level ), "Memory was not allocated for level " << level << "." );
         auto data = fmem->getPointer( level );
         return data;
      }
   }

   /// \brief Returns a pointer to the array that stores the ghost-layer data.
   ///
   /// You better know what you are doing when calling this ...
   ///
   /// \param primitiveID macro-primitive (volume primitive) to get the data of
   /// \return pointer to the array that stores the data
   ValueType* glMemory( PrimitiveID primitiveID, uint_t level, uint_t glId ) const
   {
      if ( this->storage_->hasGlobalCells() )
      {
         WALBERLA_CHECK( this->storage_->cellExistsLocally( primitiveID ),
                         "Cannot read/write DoF since macro-cell does not exists (locally)." );
         FunctionMemory< ValueType >* fmem = this->storage_->getCell( primitiveID )->getData( cellGhostLayerDataIDs_.at( glId ) );
         WALBERLA_CHECK( fmem->hasLevel( level ), "Memory was not allocated for level " << level << "." );
         ValueType* data = fmem->getPointer( level );
         return data;
      }
      else
      {
         WALBERLA_CHECK( this->storage_->faceExistsLocally( primitiveID ),
                         "Cannot read/write DoF since macro-face does not exists (locally)." );
         FunctionMemory< ValueType >* fmem = this->storage_->getFace( primitiveID )->getData( faceGhostLayerDataIDs_.at( glId ) );
         WALBERLA_CHECK( fmem->hasLevel( level ), "Memory was not allocated for level " << level << "." );
         ValueType* data = fmem->getPointer( level );
         return data;
      }
   }

   /// \brief Read access to a degree of freedom on a macro-face (for debugging / testing).
   ///
   /// This way of accessing the DoFs introduces indirections for every call. So better not use this in performance critical code.
   /// When iterating over a macro-face in performance critical code, it's better to first get the data pointer and _then_ iterate
   /// over the memory.
   ///
   /// \param primitiveID volume primitive to get the data from
   /// \param idx         index of the micro-volume
   /// \param dofID       index of the DoF on the micro-volume
   /// \param faceType    face type of the micro-volume
   /// \param level       refinement level
   /// \return value of the DoF
   ValueType
       dof( PrimitiveID primitiveID, hyteg::indexing::Index idx, uint_t dofID, facedof::FaceType faceType, uint_t level ) const
   {
      auto data = dofMemory( primitiveID, level );
      return data[indexing::index(
          idx.x(), idx.y(), faceType, dofID, numScalarsPerPrimitive_.at( primitiveID ), level, memoryLayout_ )];
   }

   /// \brief Access to a degree of freedom on a macro-face (for debugging / testing).
   ///
   /// This way of accessing the DoFs introduces indirections for every call. So better not use this in performance critical code.
   /// When iterating over a macro-face in performance critical code, it's better to first get the data pointer and _then_ iterate
   /// over the memory.
   ///
   /// \param primitiveID volume primitive to get the data from
   /// \param idx         index of the micro-volume
   /// \param dofID       index of the DoF on the micro-volume
   /// \param faceType    face type of the micro-volume
   /// \param level       refinement level
   /// \return reference to the DoF
   ValueType& dof( PrimitiveID primitiveID, hyteg::indexing::Index idx, uint_t dofID, facedof::FaceType faceType, uint_t level )
   {
      auto data = dofMemory( primitiveID, level );
      return data[indexing::index(
          idx.x(), idx.y(), faceType, dofID, numScalarsPerPrimitive_.at( primitiveID ), level, memoryLayout_ )];
   }

   /// \brief Read access to a degree of freedom on a macro-cell (for debugging / testing).
   ///
   /// This way of accessing the DoFs introduces indirections for every call. So better not use this in performance critical code.
   /// When iterating over a macro-face in performance critical code, it's better to first get the data pointer and _then_ iterate
   /// over the memory.
   ///
   /// \param primitiveID volume primitive to get the data from
   /// \param idx         index of the micro-volume
   /// \param dofID       index of the DoF on the micro-volume
   /// \param cellType    cell type of the micro-volume
   /// \param level       refinement level
   /// \return value of the DoF
   ValueType
       dof( PrimitiveID primitiveID, hyteg::indexing::Index idx, uint_t dofID, celldof::CellType cellType, uint_t level ) const
   {
      auto data = dofMemory( primitiveID, level );
      return data[indexing::index(
          idx.x(), idx.y(), idx.z(), cellType, dofID, numScalarsPerPrimitive_.at( primitiveID ), level, memoryLayout_ )];
   }

   /// \brief Access to a degree of freedom on a macro-cell (for debugging / testing).
   ///
   /// This way of accessing the DoFs introduces indirections for every call. So better not use this in performance critical code.
   /// When iterating over a macro-face in performance critical code, it's better to first get the data pointer and _then_ iterate
   /// over the memory.
   ///
   /// \param primitiveID volume primitive to get the data from
   /// \param idx         index of the micro-volume
   /// \param dofID       index of the DoF on the micro-volume
   /// \param cellType    cell type of the micro-volume
   /// \param level       refinement level
   /// \return reference to the DoF
   ValueType& dof( PrimitiveID primitiveID, hyteg::indexing::Index idx, uint_t dofID, celldof::CellType cellType, uint_t level )
   {
      auto data = dofMemory( primitiveID, level );
      return data[indexing::index(
          idx.x(), idx.y(), idx.z(), cellType, dofID, numScalarsPerPrimitive_.at( primitiveID ), level, memoryLayout_ )];
   }

   indexing::VolumeDoFMemoryLayout memoryLayout() const { return memoryLayout_; }

 private:
   using Function< VolumeDoFFunction< ValueType > >::communicators_;

   void allocateMemory();

   /// \brief Updates the map holding the number of scalars per primitive with neighborhood information.
   void communicateNumScalarsPerPrimitive();

   void addPackInfos();

   inline void deleteFunctionMemory()
   {
      if ( this->storage_->hasGlobalCells() )
      {
         this->storage_->deleteCellData( cellInnerDataID_ );
         for ( auto it : cellGhostLayerDataIDs_ )
         {
            this->storage_->deleteCellData( it.second );
         }
      }
      else
      {
         this->storage_->deleteFaceData( faceInnerDataID_ );
         for ( auto it : faceGhostLayerDataIDs_ )
         {
            this->storage_->deleteFaceData( it.second );
         }
      }
   }

   std::map< PrimitiveID, uint_t > numScalarsPerPrimitive_;

   indexing::VolumeDoFMemoryLayout memoryLayout_;

   /// Each data ID stores the inner scalar fields for all levels.
   PrimitiveDataID< FunctionMemory< ValueType >, Face > faceInnerDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Cell > cellInnerDataID_;

   /// One data ID per ghost-layer. Should be up to 3 in 2D and 4 in 3D.
   /// Maps from the local macro-edge ID (for 2D) or local macro-face ID (for 3D) to the respective ghost-layer memory.
   /// Note that there is only a single ghost-layer structure per side of the volume, even if the neighboring volume is
   /// on a different (macro-)refinement level.
   /// That means that the sizes of the ghost-layer memory might be different than the "length" of a macro-boundary on the inside.
   std::map< uint_t, PrimitiveDataID< FunctionMemory< ValueType >, Face > > faceGhostLayerDataIDs_;
   std::map< uint_t, PrimitiveDataID< FunctionMemory< ValueType >, Cell > > cellGhostLayerDataIDs_;

   /// All functions that actually allocate data and are not composites are handles to the allocated memory.
   /// This means the copy-ctor and copy-assignment only create a handle that is associated with the same memory.
   /// Deep copies must be created explicitly.
   /// To make sure that functions that are not used anymore are deleted, we need to add this reference counter to the handle.
   /// Once it drops to zero, we can deallocate the memory from the storage.
   std::shared_ptr< internal::ReferenceCounter > referenceCounter_;
};

/// \brief Given an affine coordinate (in computational space) this function locates the micro-element in the macro-face and
/// the coordinates (in reference space) local to the micro-element.
///
/// \param level            [in] refinement level
/// \param face             [in] macro-face object
/// \param coordinates      [in] the point that the local element shall include
/// \param elementIndex     [out] logical index of the micro-element
/// \param faceType         [out] face type (up- or down-face / blue or grey face)
/// \param localCoordinates [out] coordinates of the queried point translated to element local coordinates (reference space)
template < typename ValueType >
inline void getLocalElementFromCoordinates( uint_t                  level,
                                            const Face&             face,
                                            const Point2D&          coordinates,
                                            hyteg::indexing::Index& elementIndex,
                                            facedof::FaceType&      faceType,
                                            Point2D&                localCoordinates )
{
   // TODO: precompute and store trafo matrix from affine to ref space in macro-face
   //       and/or write method in macro-face that performs the trafo altogether

   // Transform absolute coordinates to macro element relative coordinates
   Matrix2r A;
   A( 0, 0 )          = ( face.getCoordinates()[1] - face.getCoordinates()[0] )[0];
   A( 0, 1 )          = ( face.getCoordinates()[2] - face.getCoordinates()[0] )[0];
   A( 1, 0 )          = ( face.getCoordinates()[1] - face.getCoordinates()[0] )[1];
   A( 1, 1 )          = ( face.getCoordinates()[2] - face.getCoordinates()[0] )[1];
   Matrix2r transform = A.inverse();

   Point2D x( { coordinates[0] - face.getCoordinates()[0][0], coordinates[1] - face.getCoordinates()[0][1] } );

   Point2D xRelMacro = transform.mul( x );

   // Determine lower-left corner index of the quad where the evaluation point lies in
   uint_t rowsize = levelinfo::num_microvertices_per_edge( level );
   real_t hInv    = walberla::real_c( rowsize - 1 );
   real_t h       = walberla::real_c( 1.0 / hInv );
   int    binX    = static_cast< int >( std::floor( xRelMacro[0] * ( rowsize - 1 ) ) );
   int    binY    = static_cast< int >( std::floor( xRelMacro[1] * ( rowsize - 1 ) ) );

   if ( binX < 0 )
   {
      binX = 0;
   }

   if ( binY < 0 )
   {
      binY = 0;
   }

   if ( binX >= static_cast< int >( rowsize - 1 ) )
   {
      binX = static_cast< int >( rowsize - 2 );
   }

   if ( binY >= static_cast< int >( rowsize - 1 ) )
   {
      binY = static_cast< int >( rowsize - 2 );
   }

   if ( binX + binY >= static_cast< int >( rowsize - 1 ) )
   {
      int binXDec = ( binX + binY - static_cast< int >( rowsize - 2 ) ) / 2;
      binXDec += ( binX + binY - static_cast< int >( rowsize - 2 ) ) % 2;
      int binYDec = ( binX + binY - static_cast< int >( rowsize - 2 ) ) / 2;
      binX -= binXDec;
      binY -= binYDec;
   }

   hyteg::indexing::Index index;
   index.x() = uint_c( binX );
   index.y() = uint_c( binY );

   WALBERLA_ASSERT_LESS( index.x(), rowsize - 1 );
   WALBERLA_ASSERT_LESS( index.y(), rowsize - 1 );
   WALBERLA_ASSERT_LESS( index.x() + index.y(), rowsize - 1, "index.x(): " << index.x() << ", index.y()" << index.y() );

   Point2D localCoords;
   localCoords[0] = xRelMacro[0] - index.x() * h;
   localCoords[1] = xRelMacro[1] - index.y() * h;
   localCoords *= hInv;

   // decide if up or down triangle
   // clamp to macro-face if the corresponding down-triangle would be out of the macro-face
   // otherwise check floating point distance
   bool upTriangle = ( index.x() + index.y() == rowsize - 2 ) || ( localCoordinates[0] + localCoordinates[1] <= 1.0 );

   elementIndex     = index;
   faceType         = upTriangle ? facedof::FaceType::GRAY : facedof::FaceType::BLUE;
   localCoordinates = localCoords;
}

/// \brief Given an affine coordinate (in computational space) this function locates the micro-element in the macro-cell and
/// the coordinates (in reference space) local to the micro-element.
///
/// \param level            [in] refinement level
/// \param cell             [in] macro-cell object
/// \param coordinates      [in] the point that the local element shall include
/// \param elementIndex     [out] logical index of the micro-element
/// \param cellType         [out] cell type
/// \param localCoordinates [out] coordinates of the queried point translated to element local coordinates (reference space)
template < typename ValueType >
inline void getLocalElementFromCoordinates( uint_t                  level,
                                            const Cell&             cell,
                                            const Point3D&          coordinates,
                                            hyteg::indexing::Index& elementIndex,
                                            celldof::CellType&      cellType,
                                            Point3D&                localCoordinates )
{
   using hyteg::indexing::Index;
   using namespace vertexdof::macrocell;

   // Assuming now that the passed coordinates are in the cell.
   // Otherwise they are clamped.

   // 1. Affine transformation to local macro-tet.

   auto xRelMacro = detail::transformToLocalTet(
       cell.getCoordinates()[0], cell.getCoordinates()[1], cell.getCoordinates()[2], cell.getCoordinates()[3], coordinates );

   // 2. Find micro-cube in macro-cell. Each micro-cube is composed of 6 cells of all 6 cell-types.

   const int    numMicroEdges = (int) levelinfo::num_microedges_per_edge( level );
   const real_t microEdgeSize = 1.0 / real_c( numMicroEdges );

   int planeX = (int) ( xRelMacro[0] / microEdgeSize );
   int planeY = (int) ( xRelMacro[1] / microEdgeSize );
   int planeZ = (int) ( xRelMacro[2] / microEdgeSize );

   // clip to prevent element outside macro-cell. std::clamp is C++17 ...
   planeX = detail::clamp( planeX, 0, numMicroEdges - 1 );
   planeY = detail::clamp( planeY, 0, numMicroEdges - 1 - planeX );
   planeZ = detail::clamp( planeZ, 0, numMicroEdges - 1 - planeX - planeY );

   // 3. In the interior of the macro-cell, the micro-cubes are contained entirely.
   //    On the boundary, some micro-tets in the micro-cube are outside of the macro-cell.
   //    Let's check that.

   // check if the point is located in a sub-cube, i.e. all cell types are possible
   const bool inFullCube = planeX + planeY + planeZ < numMicroEdges - 2;
   // check if cube is at boundary and lacks the WHITE_DOWN cell
   const bool inCutCube = planeX + planeY + planeZ == numMicroEdges - 2;
   // if both are false, the point is located in a white cell near the macro-edges / macro-vertices

   // 4. Now we got through all <= 6 micro-tets in the micro-cube and check if the point is contained.
   //    We can perform some shortcuts, in general, however, we perform the coordinate transformation
   //    for all cells.
   //    If no micro-tet can be chosen directly (due to floating-point math errors) we choose the "closest"
   //    micro-tet by minimal summed distance to all faces.
   //
   //    Possible optimization: the transformation to the 6 local tets can be optimized by transforming to the
   //                           local micro-cube space first, and then apply the precomputed transforms to the point.
   //                           That could save up to 5 3x3 matrix solves(!)

   if ( ( !inFullCube && !inCutCube ) )
   {
      cellType = celldof::CellType::WHITE_UP;
   }
   else
   {
      std::vector< celldof::CellType > possibleCellTypes = { celldof::CellType::WHITE_UP,
                                                             celldof::CellType::BLUE_UP,
                                                             celldof::CellType::GREEN_UP,
                                                             celldof::CellType::BLUE_DOWN,
                                                             celldof::CellType::GREEN_DOWN };
      if ( inFullCube )
      {
         possibleCellTypes.push_back( celldof::CellType::WHITE_DOWN );
      }

      real_t maxDistSum = std::numeric_limits< real_t >::max();

      for ( auto ct : possibleCellTypes )
      {
         auto mci = celldof::macrocell::getMicroVerticesFromMicroCell( Index( planeX, planeY, planeZ ), ct );
         auto mt0 = coordinateFromIndex( level, cell, mci[0] );
         auto mt1 = coordinateFromIndex( level, cell, mci[1] );
         auto mt2 = coordinateFromIndex( level, cell, mci[2] );
         auto mt3 = coordinateFromIndex( level, cell, mci[3] );

         auto xl = detail::transformToLocalTet( mt0, mt1, mt2, mt3, coordinates );
         auto s  = xl[0] + xl[1] + xl[2];

         Point4D rel( { xl[0], xl[1], xl[2], s } );

         real_t distSum  = 0;
         bool   contains = true;
         for ( uint_t i = 0; i < 4; i++ )
         {
            if ( rel[i] < 0 )
            {
               distSum += std::abs( rel[i] );
               contains = false;
            }
            else if ( rel[i] > 1 )
            {
               distSum += std::abs( rel[i] - 1 );
               contains = false;
            }
         }

         if ( contains )
         {
            cellType = ct;
            break;
         }

         if ( distSum < maxDistSum )
         {
            cellType   = ct;
            maxDistSum = distSum;
         }
      }
   }

   auto microCellIndices = celldof::macrocell::getMicroVerticesFromMicroCell( Index( planeX, planeY, planeZ ), cellType );

   auto microTet0 = vertexdof::macrocell::coordinateFromIndex( level, cell, microCellIndices[0] );
   auto microTet1 = vertexdof::macrocell::coordinateFromIndex( level, cell, microCellIndices[1] );
   auto microTet2 = vertexdof::macrocell::coordinateFromIndex( level, cell, microCellIndices[2] );
   auto microTet3 = vertexdof::macrocell::coordinateFromIndex( level, cell, microCellIndices[3] );

   localCoordinates =
       vertexdof::macrocell::detail::transformToLocalTet( microTet0, microTet1, microTet2, microTet3, coordinates );

   elementIndex.x() = planeX;
   elementIndex.y() = planeY;
   elementIndex.z() = planeZ;

   const uint_t numMicroVertices = levelinfo::num_microvertices_per_edge( level );
   WALBERLA_ASSERT_LESS( microCellIndices[0].x() + microCellIndices[0].y() + microCellIndices[0].z(), numMicroVertices );
   WALBERLA_ASSERT_LESS( microCellIndices[1].x() + microCellIndices[1].y() + microCellIndices[1].z(), numMicroVertices );
   WALBERLA_ASSERT_LESS( microCellIndices[2].x() + microCellIndices[2].y() + microCellIndices[2].z(), numMicroVertices );
   WALBERLA_ASSERT_LESS( microCellIndices[3].x() + microCellIndices[3].y() + microCellIndices[3].z(), numMicroVertices );

   WALBERLA_DEBUG_SECTION()
   {
      auto xLocal = localCoordinates;

      auto sum = xLocal[0] + xLocal[1] + xLocal[2];

      if ( xLocal[0] < -1e-8 || xLocal[0] > 1.0 + 1e-8 || xLocal[1] < -1e-8 || xLocal[1] > 1.0 + 1e-8 || xLocal[2] < -1e-8 ||
           xLocal[2] > 1.0 + 1e-8 || sum < -1e-8 || sum > 1.0 + 1e-8 )
      {
         WALBERLA_LOG_DEVEL( "Bad cell choice:" )
         WALBERLA_LOG_DEVEL( "local " << xLocal << ", global " << coordinates
                                      << ", local sum: " << walberla::format( "%10.5e", sum ) << ", micro-cell type"
                                      << celldof::CellTypeToStr.at( cellType ) );
         for ( auto ct : celldof::allCellTypes )
         {
            auto mci = celldof::macrocell::getMicroVerticesFromMicroCell( Index( planeX, planeY, planeZ ), ct );
            auto mt0 = coordinateFromIndex( level, cell, mci[0] );
            auto mt1 = coordinateFromIndex( level, cell, mci[1] );
            auto mt2 = coordinateFromIndex( level, cell, mci[2] );
            auto mt3 = coordinateFromIndex( level, cell, mci[3] );

            auto xl = detail::transformToLocalTet( mt0, mt1, mt2, mt3, coordinates );
            auto s  = xl[0] + xl[1] + xl[2];
            WALBERLA_LOG_DEVEL( "Cell type: " << celldof::CellTypeToStr.at( ct ) );
            WALBERLA_LOG_DEVEL( "x local:   " << xl << ", local sum: " << walberla::format( "%10.5e", s ) );
         }

         WALBERLA_LOG_DEVEL( "Breakpoint here ..." );
      }
   }
}

} // namespace volumedofspace
} // namespace hyteg
