/*
 * Copyright (c) 2025 Marcus Mohr.
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

#include <algorithm>
#include <vector>

#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivedata/PrimitiveDataID.hpp"
#include "hyteg/primitives/Primitive.hpp"
#include "hyteg/primitives/PrimitiveID.hpp"
#include "hyteg/primitives/all.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/types/Concepts.hpp"

namespace hyteg {

class SubStorage
{
 public:
   /// Function to create a subset of a SetupPrimitiveStorage from an existing SetupPrimitiveStorage using an oracle for filtering
   template < typename T >
      requires requires( T                                      oracle,
                         const std::shared_ptr< const Vertex >& vertex,
                         const std::shared_ptr< const Edge >&   edge,
                         const std::shared_ptr< const Face >&   face,
                         const std::shared_ptr< const Cell >&   cell ) {
         {
            oracle( vertex )
         } -> std::same_as< bool >;

         {
            oracle( edge )
         } -> std::same_as< bool >;

         {
            oracle( face )
         } -> std::same_as< bool >;

         {
            oracle( cell )
         } -> std::same_as< bool >;
      }
   static std::shared_ptr< SetupPrimitiveStorage > extractSubsetFromSetupStorage( const SetupPrimitiveStorage& setupStore,
                                                                                  T&&                          oracle )
   {
      // STEP 1:
      //
      // We loop over all primitives in the SetupPrimitiveStorage to figure out, if they belong to
      // the desired subset, or not.
      std::vector< PrimitiveID > allSubIDs;
      bool                       subsetHasCells = false;

      auto getSubsetPids = [&allSubIDs, &oracle]( const auto& map ) {
         for ( const auto& item : map )
         {
            if ( oracle( item.second ) )
            {
               allSubIDs.push_back( item.first );
            }
         }
      };

      getSubsetPids( setupStore.getCells() );
      getSubsetPids( setupStore.getFaces() );
      getSubsetPids( setupStore.getEdges() );
      getSubsetPids( setupStore.getVertices() );

      // STEP 2:
      //
      // Create clones of all primitives belonging to the subset and adapt their neighborhood information
      SetupPrimitiveStorage::CellMap   subCellMap;
      SetupPrimitiveStorage::FaceMap   subFaceMap;
      SetupPrimitiveStorage::EdgeMap   subEdgeMap;
      SetupPrimitiveStorage::VertexMap subVertexMap;

      // return true, if testPid is _not_ in list of subPids
      auto pidNotInSubset = [&allSubIDs]( const PrimitiveID& testPid ) -> bool {
         return ( find( allSubIDs.begin(), allSubIDs.end(), testPid ) == allSubIDs.end() );
      };

      // return true, if testPid is _not_ in list of subPids
      auto pairPidNotInSubset = [&allSubIDs]( const std::pair< uint_t, PrimitiveID >& testPair ) -> bool {
         return ( find( allSubIDs.begin(), allSubIDs.end(), testPair.second ) == allSubIDs.end() );
      };

      for ( const auto& item : setupStore.getVertices() )
      {
         PrimitiveID pid = item.first;
         if ( std::find( allSubIDs.begin(), allSubIDs.end(), pid ) != allSubIDs.end() )
         {
            // create a clone
            auto newVertex = std::make_shared< Vertex >( *( item.second ) );

            // adapt clone
            newVertex->meshBoundaryFlag_ = 0;
            std::erase_if( newVertex->neighborCells_, pidNotInSubset );
            std::erase_if( newVertex->neighborFaces_, pidNotInSubset );
            std::erase_if( newVertex->neighborEdges_, pidNotInSubset );
            std::erase_if( newVertex->neighborVertices_, pidNotInSubset );

            // insert clone into map
            subVertexMap[pid] = newVertex;
         }
      }

      for ( const auto& item : setupStore.getEdges() )
      {
         PrimitiveID pid = item.first;
         if ( std::find( allSubIDs.begin(), allSubIDs.end(), pid ) != allSubIDs.end() )
         {
            // create a clone
            auto newEdge = std::make_shared< Edge >( *( item.second ) );

            // adapt clone
            newEdge->meshBoundaryFlag_ = 0;
            std::erase_if( newEdge->neighborCells_, pidNotInSubset );
            std::erase_if( newEdge->neighborFaces_, pidNotInSubset );
            std::erase_if( newEdge->neighborEdges_, pidNotInSubset );
            std::erase_if( newEdge->neighborVertices_, pidNotInSubset );

            // insert clone into map
            subEdgeMap[pid] = newEdge;
         }
      }

      for ( const auto& item : setupStore.getFaces() )
      {
         PrimitiveID pid = item.first;
         if ( std::find( allSubIDs.begin(), allSubIDs.end(), pid ) != allSubIDs.end() )
         {
            // create a clone
            auto newFace = std::make_shared< Face >( *( item.second ) );

            // adapt clone
            newFace->meshBoundaryFlag_ = 0;
            std::erase_if( newFace->neighborCells_, pidNotInSubset );
            std::erase_if( newFace->neighborFaces_, pidNotInSubset );
            std::erase_if( newFace->neighborEdges_, pidNotInSubset );
            std::erase_if( newFace->neighborVertices_, pidNotInSubset );

            // correct the extra information faces carry
            std::erase_if( newFace->indirectNeighborFaceIDsOverVertices_, pidNotInSubset );
            std::erase_if( newFace->indirectNeighborFaceIDsOverEdges_, pairPidNotInSubset );

            // insert clone into map
            subFaceMap[pid] = newFace;
         }
      }

      for ( const auto& item : setupStore.getCells() )
      {
         PrimitiveID pid = item.first;
         if ( std::find( allSubIDs.begin(), allSubIDs.end(), pid ) != allSubIDs.end() )
         {
            // create a clone
            auto newCell = std::make_shared< Cell >( *( item.second ) );

            // adapt clone
            newCell->meshBoundaryFlag_ = 0;
            std::erase_if( newCell->neighborCells_, pidNotInSubset );
            std::erase_if( newCell->neighborFaces_, pidNotInSubset );
            std::erase_if( newCell->neighborEdges_, pidNotInSubset );
            std::erase_if( newCell->neighborVertices_, pidNotInSubset );

            // correct the extra information faces carry
            std::erase_if( newCell->indirectNeighborCellIDsOverVertices_, pidNotInSubset );
            std::erase_if( newCell->indirectNeighborCellIDsOverFaces_, pairPidNotInSubset );

            // insert clone into map
            subCellMap[pid] = newCell;
         }
      }

      // STEP 3:
      //
      // We want the cloned sub-primitives to be on the same MPI process as the original primitives
      auto primitiveIDToTargetRankMap = std::make_shared< std::map< PrimitiveID, uint_t > >();
      for ( const auto& pid : allSubIDs )
      {
         primitiveIDToTargetRankMap->insert( { pid, setupStore.getTargetRank( pid ) } );
      }

      // STEP 4:
      //
      // Create the extrated SetupPrimitiveStorage object
      auto subStore = std::make_shared< SetupPrimitiveStorage >( subVertexMap,
                                                                 subEdgeMap,
                                                                 subFaceMap,
                                                                 subCellMap,
                                                                 uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ),
                                                                 primitiveIDToTargetRankMap );

      // std::shared_ptr< SetupPrimitiveStorage > subStore;

      return subStore;
   }

   /// An assign function that works for P1Functions (potentially) defined on different PrimitiveStorage objects
   ///
   /// The function works in the same fashion as the standard assign() of FEFunctions. However, there is a difference,
   /// it will only perform assignment for the primitives inside the targetStorage. Consequently the PrimitiveStorage
   /// of the dstFunction and all srcFunctions must either be identical to the targetStorage or be a superset.
   template < concepts::value_type ValueType >
   static void assignAcrossStorages( const std::shared_ptr< PrimitiveStorage >&                                    targetStorage,
                                     const std::vector< ValueType >&                                               scalars,
                                     const P1Function< ValueType >&                                                dstFunction,
                                     const std::vector< std::reference_wrapper< const P1Function< ValueType > > >& srcFunctions,
                                     uint_t                                                                        level,
                                     DoFType                                                                       flag = All )
   {
      WALBERLA_CHECK_EQUAL( scalars.size(), srcFunctions.size() )

      // Collect all source IDs in a vector
      std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
      std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
      std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;
      std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >   srcCellIDs;

      for ( const P1Function< ValueType >& function : srcFunctions )
      {
         srcVertexIDs.push_back( function.getVertexDataID() );
         srcEdgeIDs.push_back( function.getEdgeDataID() );
         srcFaceIDs.push_back( function.getFaceDataID() );
         srcCellIDs.push_back( function.getCellDataID() );
      }

      // IDs for access to destination buffers
      PrimitiveDataID< FunctionMemory< ValueType >, Vertex > dstVertexID = dstFunction.getVertexDataID();
      PrimitiveDataID< FunctionMemory< ValueType >, Edge >   dstEdgeID   = dstFunction.getEdgeDataID();
      PrimitiveDataID< FunctionMemory< ValueType >, Face >   dstFaceID   = dstFunction.getFaceDataID();
      PrimitiveDataID< FunctionMemory< ValueType >, Cell >   dstCellID   = dstFunction.getCellDataID();

      // --------------------
      //  Assign on Vertices
      // --------------------
      std::vector< std::shared_ptr< Vertex > > srcVertices;
      auto                                     numVertices = targetStorage->getVertices().size();
      for ( auto& item : targetStorage->getVertices() )
      {
         PrimitiveID pid = item.first;

         // need to be careful here, we must ensure that we take the vertex object from the
         // PrimitiveStorage that belongs to the dstFunction! This might be different from
         // item.second!
         std::shared_ptr< Vertex > dstVertex = dstFunction.getStorage()->getSharedPointerToVertex( pid );
         WALBERLA_ASSERT_NOT_NULLPTR( dstVertex );

         // assemble list of vertices from which we will obtain the rhs data
         srcVertices.clear();
         for ( const auto& func : srcFunctions )
         {
            srcVertices.push_back( func.get().getStorage()->getSharedPointerToVertex( pid ) );
            WALBERLA_ASSERT_NOT_NULLPTR( func.get().getStorage()->getSharedPointerToVertex( pid ) );
         }

         if ( testFlag( dstFunction.getBoundaryCondition().getBoundaryType( dstVertex->getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macrovertex::assignAcrossStorages< ValueType >(
                srcVertices, dstVertex, scalars, srcVertexIDs, dstVertexID, level );
         }
      }

      // -----------------
      //  Assign on Edges
      // -----------------
      std::vector< std::shared_ptr< Edge > > srcEdges;
      for ( auto& item : targetStorage->getEdges() )
      {
         PrimitiveID pid = item.first;

         // need to be careful here, we must ensure that we take the edge object from the
         // PrimitiveStorage that belongs to the dstFunction! This might be different from
         // item.second!
         std::shared_ptr< Edge > dstEdge = dstFunction.getStorage()->getSharedPointerToEdge( pid );
         WALBERLA_ASSERT_NOT_NULLPTR( dstEdge );

         // assemble list of edge from which we will obtain the rhs data
         srcEdges.clear();
         for ( const auto& func : srcFunctions )
         {
            srcEdges.push_back( func.get().getStorage()->getSharedPointerToEdge( pid ) );
            WALBERLA_ASSERT_NOT_NULLPTR( func.get().getStorage()->getSharedPointerToEdge( pid ) );
         }

         if ( testFlag( dstFunction.getBoundaryCondition().getBoundaryType( dstEdge->getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macroedge::assignAcrossStorages< ValueType >( srcEdges, dstEdge, scalars, srcEdgeIDs, dstEdgeID, level );
         }
      }

      // -----------------
      //  Assign on Faces
      // -----------------
      std::vector< std::shared_ptr< Face > > srcFaces;
      for ( auto& item : targetStorage->getFaces() )
      {
         PrimitiveID pid = item.first;

         // need to be careful here, we must ensure that we take the face object from the
         // PrimitiveStorage that belongs to the dstFunction! This might be different from
         // item.second!
         std::shared_ptr< Face > dstFace = dstFunction.getStorage()->getSharedPointerToFace( pid );
         WALBERLA_ASSERT_NOT_NULLPTR( dstFace );

         // assemble list of edge from which we will obtain the rhs data
         srcFaces.clear();
         for ( const auto& func : srcFunctions )
         {
            srcFaces.push_back( func.get().getStorage()->getSharedPointerToFace( pid ) );
            WALBERLA_ASSERT_NOT_NULLPTR( func.get().getStorage()->getSharedPointerToFace( pid ) );
         }

         if ( testFlag( dstFunction.getBoundaryCondition().getBoundaryType( dstFace->getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macroface::assignAcrossStorages< ValueType >( srcFaces, dstFace, scalars, srcFaceIDs, dstFaceID, level );
         }
      }

      // -----------------
      //  Assign on Cells
      // -----------------
      std::vector< std::shared_ptr< Cell > > srcCells;
      for ( auto& item : targetStorage->getCells() )
      {
         PrimitiveID pid = item.first;

         // need to be careful here, we must ensure that we take the cell object from the
         // PrimitiveStorage that belongs to the dstFunction! This might be different from
         // item.second!
         std::shared_ptr< Cell > dstCell = dstFunction.getStorage()->getSharedPointerToCell( pid );
         WALBERLA_ASSERT_NOT_NULLPTR( dstCell );

         // assemble list of edge from which we will obtain the rhs data
         srcCells.clear();
         for ( const auto& func : srcFunctions )
         {
            srcCells.push_back( func.get().getStorage()->getSharedPointerToCell( pid ) );
            WALBERLA_ASSERT_NOT_NULLPTR( func.get().getStorage()->getSharedPointerToCell( pid ) );
         }

         if ( testFlag( dstFunction.getBoundaryCondition().getBoundaryType( dstCell->getMeshBoundaryFlag() ), flag ) )
         {
            vertexdof::macrocell::assignAcrossStorages< ValueType >( srcCells, dstCell, scalars, srcCellIDs, dstCellID, level );
         }
      }
   }

   /// An assign function that works for EdgeDoFFunctions (potentially) defined on different PrimitiveStorage objects
   ///
   /// The function works in the same fashion as the standard assign() of FEFunctions. However, there is a difference,
   /// it will only perform assignment for the primitives inside the targetStorage. Consequently the PrimitiveStorage
   /// of the dstFunction and all srcFunctions must either be identical to the targetStorage or be a superset.
   template < concepts::value_type ValueType >
   static void
       assignAcrossStorages( const std::shared_ptr< PrimitiveStorage >&                                         targetStorage,
                             const std::vector< ValueType >&                                                    scalars,
                             const EdgeDoFFunction< ValueType >&                                                dstFunction,
                             const std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >& srcFunctions,
                             uint_t                                                                             level,
                             DoFType                                                                            flag = All )
   {
      WALBERLA_CHECK_EQUAL( scalars.size(), srcFunctions.size() )

      // Collect all source IDs in a vector
      std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > > srcEdgeIDs;
      std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > srcFaceIDs;
      std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > srcCellIDs;

      for ( const EdgeDoFFunction< ValueType >& function : srcFunctions )
      {
         srcEdgeIDs.push_back( function.getEdgeDataID() );
         srcFaceIDs.push_back( function.getFaceDataID() );
         srcCellIDs.push_back( function.getCellDataID() );
      }

      // IDs for access to destination buffers
      PrimitiveDataID< FunctionMemory< ValueType >, Edge > dstEdgeID = dstFunction.getEdgeDataID();
      PrimitiveDataID< FunctionMemory< ValueType >, Face > dstFaceID = dstFunction.getFaceDataID();
      PrimitiveDataID< FunctionMemory< ValueType >, Cell > dstCellID = dstFunction.getCellDataID();

      // -----------------
      //  Assign on Edges
      // -----------------
      std::vector< std::shared_ptr< Edge > > srcEdges;
      for ( auto& item : targetStorage->getEdges() )
      {
         PrimitiveID pid = item.first;

         // need to be careful here, we must ensure that we take the edge object from the
         // PrimitiveStorage that belongs to the dstFunction! This might be different from
         // item.second!
         std::shared_ptr< Edge > dstEdge = dstFunction.getStorage()->getSharedPointerToEdge( pid );
         WALBERLA_ASSERT_NOT_NULLPTR( dstEdge );

         // assemble list of edge from which we will obtain the rhs data
         srcEdges.clear();
         for ( const auto& func : srcFunctions )
         {
            srcEdges.push_back( func.get().getStorage()->getSharedPointerToEdge( pid ) );
            WALBERLA_ASSERT_NOT_NULLPTR( func.get().getStorage()->getSharedPointerToEdge( pid ) );
         }

         if ( testFlag( dstFunction.getBoundaryCondition().getBoundaryType( dstEdge->getMeshBoundaryFlag() ), flag ) )
         {
            edgedof::macroedge::assignAcrossStorages< ValueType >( srcEdges, dstEdge, scalars, srcEdgeIDs, dstEdgeID, level );
         }
      }

      // -----------------
      //  Assign on Faces
      // -----------------
      std::vector< std::shared_ptr< Face > > srcFaces;
      for ( auto& item : targetStorage->getFaces() )
      {
         PrimitiveID pid = item.first;

         // need to be careful here, we must ensure that we take the face object from the
         // PrimitiveStorage that belongs to the dstFunction! This might be different from
         // item.second!
         std::shared_ptr< Face > dstFace = dstFunction.getStorage()->getSharedPointerToFace( pid );
         WALBERLA_ASSERT_NOT_NULLPTR( dstFace );

         // assemble list of edge from which we will obtain the rhs data
         srcFaces.clear();
         for ( const auto& func : srcFunctions )
         {
            srcFaces.push_back( func.get().getStorage()->getSharedPointerToFace( pid ) );
            WALBERLA_ASSERT_NOT_NULLPTR( func.get().getStorage()->getSharedPointerToFace( pid ) );
         }

         if ( testFlag( dstFunction.getBoundaryCondition().getBoundaryType( dstFace->getMeshBoundaryFlag() ), flag ) )
         {
            edgedof::macroface::assignAcrossStorages< ValueType >( srcFaces, dstFace, scalars, srcFaceIDs, dstFaceID, level );
         }
      }

      // -----------------
      //  Assign on Cells
      // -----------------
      std::vector< std::shared_ptr< Cell > > srcCells;
      for ( auto& item : targetStorage->getCells() )
      {
         PrimitiveID pid = item.first;

         // need to be careful here, we must ensure that we take the cell object from the
         // PrimitiveStorage that belongs to the dstFunction! This might be different from
         // item.second!
         std::shared_ptr< Cell > dstCell = dstFunction.getStorage()->getSharedPointerToCell( pid );
         WALBERLA_ASSERT_NOT_NULLPTR( dstCell );

         // assemble list of edge from which we will obtain the rhs data
         srcCells.clear();
         for ( const auto& func : srcFunctions )
         {
            srcCells.push_back( func.get().getStorage()->getSharedPointerToCell( pid ) );
            WALBERLA_ASSERT_NOT_NULLPTR( func.get().getStorage()->getSharedPointerToCell( pid ) );
         }

         if ( testFlag( dstFunction.getBoundaryCondition().getBoundaryType( dstCell->getMeshBoundaryFlag() ), flag ) )
         {
            edgedof::macrocell::assignAcrossStorages< ValueType >( srcCells, dstCell, scalars, srcCellIDs, dstCellID, level );
         }
      }
   }

   /// An assign function that works for P2Functions (potentially) defined on different PrimitiveStorage objects
   ///
   /// The function works in the same fashion as the standard assign() of FEFunctions. However, there is a difference,
   /// it will only perform assignment for the primitives inside the targetStorage. Consequently the PrimitiveStorage
   /// of the dstFunction and all srcFunctions must either be identical to the targetStorage or be a superset.
   template < concepts::value_type ValueType >
   static void assignAcrossStorages( const std::shared_ptr< PrimitiveStorage >&                                    targetStorage,
                                     const std::vector< ValueType >&                                               scalars,
                                     const P2Function< ValueType >&                                                dstFunction,
                                     const std::vector< std::reference_wrapper< const P2Function< ValueType > > >& srcFunctions,
                                     uint_t                                                                        level,
                                     DoFType                                                                       flag = All )
   {
      WALBERLA_ASSERT_EQUAL( scalars.size(), srcFunctions.size() )

      std::vector< std::reference_wrapper< const vertexdof::VertexDoFFunction< ValueType > > > srcVertexDoFFunctions;
      std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >              srcEdgeDoFFunctions;

      for ( const P2Function< ValueType >& function : srcFunctions )
      {
         srcVertexDoFFunctions.push_back( function.getVertexDoFFunction() );
         srcEdgeDoFFunctions.push_back( function.getEdgeDoFFunction() );
      }

      assignAcrossStorages( targetStorage, scalars, dstFunction.getVertexDoFFunction(), srcVertexDoFFunctions, level, flag );
      assignAcrossStorages( targetStorage, scalars, dstFunction.getEdgeDoFFunction(), srcEdgeDoFFunctions, level, flag );
   }
};

} // namespace hyteg
