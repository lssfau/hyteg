/*
 * Copyright (c) 2023 Marcus Mohr.
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

#include "hyteg/communication/Syncing.hpp"

#include "core/Environment.h"
#include "core/debug/CheckFunctions.h"
#include "core/debug/TestSubsystem.h"

#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"

using namespace hyteg;

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

const uint_t vrtxDoFValue{ 1 };
const uint_t edgeDoFValue{ 2 };
const uint_t faceDoFValue{ 3 };

template < typename value_t >
void prepareDataValues( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level, P2Function< value_t >& func )
{
   func.setToZero( level );

   const auto& vFunc = func.getVertexDoFFunction();
   const auto& eFunc = func.getEdgeDoFFunction();

   const value_t vertexConstant = static_cast< value_t >( vrtxDoFValue );
   const value_t edgeConstant   = static_cast< value_t >( edgeDoFValue );
   const value_t faceConstant   = static_cast< value_t >( faceDoFValue );

   // --------------------------------------------------
   //  set DoF values (not belonging to halos) for face
   // --------------------------------------------------
   std::vector< PrimitiveID > faceIDs = storage->getFaceIDs();
   if ( faceIDs.size() != 0 )
   {
      WALBERLA_CHECK_EQUAL( faceIDs.size(), 1u );
      Face* facePtr = dynamic_cast< Face* >( storage->getPrimitive( faceIDs[0] ) );

      value_t* faceData = facePtr->getData( vFunc.getFaceDataID() )->getPointer( level );
      for ( const auto& it : vertexdof::macroface::Iterator( level, 1 ) )
      {
         const uint_t idx = vertexdof::macroface::indexFromVertex( level, it.x(), it.y(), stencilDirection::VERTEX_C );
         faceData[idx]    = faceConstant;
      }

      faceData = facePtr->getData( eFunc.getFaceDataID() )->getPointer( level );

      for ( const auto& it : edgedof::macroface::Iterator( level, 0 ) )
      {
         // Do not update horizontal DoFs at bottom
         if ( it.y() != 0 )
         {
            faceData[edgedof::macroface::horizontalIndex( level, it.x(), it.y() )] = faceConstant;
         }

         // Do not update vertical DoFs at left border
         if ( it.x() != 0 )
         {
            faceData[edgedof::macroface::verticalIndex( level, it.x(), it.y() )] = faceConstant;
         }

         // Do not update diagonal DoFs at diagonal border
         if ( it.x() + it.y() != static_cast< idx_t >( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
         {
            faceData[edgedof::macroface::diagonalIndex( level, it.x(), it.y() )] = faceConstant;
         }
      }
   }

   // ---------------------------------------------------
   //  set DoF values (not belonging to halos) for edges
   // ---------------------------------------------------
   std::vector< PrimitiveID > edgeIDs = storage->getEdgeIDs();
   if ( edgeIDs.size() != 0 )
   {
      WALBERLA_CHECK_LESS_EQUAL( edgeIDs.size(), 3u );
      for ( const auto& iter : storage->getEdges() )
      {
         const Edge& edge  = *iter.second;
         value_t*    vDoFs = edge.getData( vFunc.getEdgeDataID() )->getPointer( level );
         value_t*    eDoFs = edge.getData( eFunc.getEdgeDataID() )->getPointer( level );

         for ( const auto& it : vertexdof::macroedge::Iterator( level, 1 ) )
         {
            vDoFs[vertexdof::macroedge::indexFromVertex( level, it.x(), stencilDirection::VERTEX_C )] = edgeConstant;
         }

         for ( const auto& it : edgedof::macroedge::Iterator( level ) )
         {
            eDoFs[edgedof::macroedge::indexFromHorizontalEdge( level, it.x(), stencilDirection::EDGE_HO_C )] = edgeConstant;
         }
      }
   }

   // ------------------------------------------------------
   //  set DoF values (not belonging to halos) for vertices
   // ------------------------------------------------------
   std::vector< PrimitiveID > vrtxIDs = storage->getVertexIDs();
   if ( vrtxIDs.size() != 0 )
   {
      WALBERLA_CHECK_LESS_EQUAL( vrtxIDs.size(), 3u );
      for ( const auto& it : storage->getVertices() )
      {
         const Vertex& vertex                                              = *it.second;
         vertex.getData( vFunc.getVertexDataID() )->getPointer( level )[0] = vertexConstant;
      }
   }
}

// return the minimal and maximal DoF value over the three different kinds of primitives in 2D;
// hereby the halos of the primitives are included
template < typename value_t >
auto getPrimitiveExtrema( const std::shared_ptr< PrimitiveStorage >& storage, uint_t level, P2Function< value_t >& func )
{
   value_t minVertex = std::numeric_limits< value_t >::max();
   value_t minEdge   = std::numeric_limits< value_t >::max();
   value_t minFace   = std::numeric_limits< value_t >::max();

   value_t maxVertex = std::numeric_limits< value_t >::min();
   value_t maxEdge   = std::numeric_limits< value_t >::min();
   value_t maxFace   = std::numeric_limits< value_t >::min();

   const auto& vFunc = func.getVertexDoFFunction();
   const auto& eFunc = func.getEdgeDoFFunction();

   for ( const auto& it : storage->getVertices() )
   {
      const Vertex& vertex = *it.second;
      for ( const auto& value : vertex.getData( vFunc.getVertexDataID() )->getVector( level ) )
      {
         minVertex = value < minVertex ? value : minVertex;
         maxVertex = value > maxVertex ? value : maxVertex;
      }
      for ( const auto& value : vertex.getData( eFunc.getVertexDataID() )->getVector( level ) )
      {
         minVertex = value < minVertex ? value : minVertex;
         maxVertex = value > maxVertex ? value : maxVertex;
      }
   }

   for ( const auto& it : storage->getEdges() )
   {
      const Edge& edge = *it.second;
      for ( const auto& value : edge.getData( vFunc.getEdgeDataID() )->getVector( level ) )
      {
         minEdge = value < minEdge ? value : minEdge;
         maxEdge = value > maxEdge ? value : maxEdge;
      }
      for ( const auto& value : edge.getData( eFunc.getEdgeDataID() )->getVector( level ) )
      {
         minEdge = value < minEdge ? value : minEdge;
         maxEdge = value > maxEdge ? value : maxEdge;
      }
   }

   for ( const auto& it : storage->getFaces() )
   {
      const Face& face = *it.second;
      for ( const auto& value : face.getData( vFunc.getFaceDataID() )->getVector( level ) )
      {
         minFace = value < minFace ? value : minFace;
         maxFace = value > maxFace ? value : maxFace;
      }
      for ( const auto& value : face.getData( eFunc.getFaceDataID() )->getVector( level ) )
      {
         minFace = value < minFace ? value : minFace;
         maxFace = value > maxFace ? value : maxFace;
      }
   }

   value_t minVertexGlobal = -walberla::mpi::allReduce( -minVertex, walberla::mpi::MAX );
   value_t minEdgeGlobal   = -walberla::mpi::allReduce( -minEdge, walberla::mpi::MAX );
   value_t minFaceGlobal   = -walberla::mpi::allReduce( -minFace, walberla::mpi::MAX );

   value_t maxVertexGlobal = walberla::mpi::allReduce( maxVertex, walberla::mpi::MAX );
   value_t maxEdgeGlobal   = walberla::mpi::allReduce( maxEdge, walberla::mpi::MAX );
   value_t maxFaceGlobal   = walberla::mpi::allReduce( maxFace, walberla::mpi::MAX );

   return std::make_tuple( minVertexGlobal, minEdgeGlobal, minFaceGlobal, maxVertexGlobal, maxEdgeGlobal, maxFaceGlobal );
}

template < typename value_t >
void runTest( const std::shared_ptr< PrimitiveStorage >& storage )
{
   const uint_t level = 3;

   P2Function< value_t > func( "P2 function", storage, level, level );

   // check values before any synchronisation
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Testing without synchronisation" );
      prepareDataValues( storage, level, func );
      auto [minVrtx, minEdge, minFace, maxVrtx, maxEdge, maxFace] = getPrimitiveExtrema( storage, level, func );

      WALBERLA_CHECK_EQUAL( minVrtx, static_cast< value_t >( 0 ) );
      WALBERLA_CHECK_EQUAL( minEdge, static_cast< value_t >( 0 ) );
      WALBERLA_CHECK_EQUAL( minFace, static_cast< value_t >( 0 ) );

      WALBERLA_CHECK_EQUAL( maxVrtx, static_cast< value_t >( vrtxDoFValue ) );
      WALBERLA_CHECK_EQUAL( maxEdge, static_cast< value_t >( edgeDoFValue ) );
      WALBERLA_CHECK_EQUAL( maxFace, static_cast< value_t >( faceDoFValue ) );
   }

   // check values after LOW2HIGH synchronisation
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Testing with LOW2HIGH synchronisation" );
      prepareDataValues( storage, level, func );
      communication::syncFunctionBetweenPrimitives( func, level, communication::syncDirection_t::LOW2HIGH );
      auto [minVrtx, minEdge, minFace, maxVrtx, maxEdge, maxFace] = getPrimitiveExtrema( storage, level, func );

      WALBERLA_CHECK_EQUAL( minVrtx, static_cast< value_t >( 0 ) );
      WALBERLA_CHECK_EQUAL( minEdge, static_cast< value_t >( 0 ) );
      WALBERLA_CHECK_EQUAL( minFace, static_cast< value_t >( vrtxDoFValue ) );

      WALBERLA_CHECK_EQUAL( maxVrtx, static_cast< value_t >( vrtxDoFValue ) );
      WALBERLA_CHECK_EQUAL( maxEdge, static_cast< value_t >( edgeDoFValue ) );
      WALBERLA_CHECK_EQUAL( maxFace, static_cast< value_t >( faceDoFValue ) );
   }

   // check values after BIDIRECTIONAL synchronisation
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Testing with BIDIRECTIONAL synchronisation" );
      prepareDataValues( storage, level, func );
      communication::syncFunctionBetweenPrimitives( func, level, communication::syncDirection_t::BIDIRECTIONAL );
      auto [minVrtx, minEdge, minFace, maxVrtx, maxEdge, maxFace] = getPrimitiveExtrema( storage, level, func );

      WALBERLA_CHECK_EQUAL( minVrtx, static_cast< value_t >( vrtxDoFValue ) );
      WALBERLA_CHECK_EQUAL( minEdge, static_cast< value_t >( vrtxDoFValue ) );
      WALBERLA_CHECK_EQUAL( minFace, static_cast< value_t >( vrtxDoFValue ) );

      WALBERLA_CHECK_EQUAL( maxVrtx, static_cast< value_t >( faceDoFValue ) );
      WALBERLA_CHECK_EQUAL( maxEdge, static_cast< value_t >( faceDoFValue ) );
      WALBERLA_CHECK_EQUAL( maxFace, static_cast< value_t >( faceDoFValue ) );
   }
}

int main( int argc, char* argv[] )
{
   walberla::debug::enterTestMode();

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   MeshInfo                            mesh = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "tri_1el.msh" ) );
   SetupPrimitiveStorage               setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   WALBERLA_LOG_INFO_ON_ROOT( "*** Value Type = double ***" );
   runTest< double >( storage );

   WALBERLA_LOG_INFO_ON_ROOT( "*** Value Type = int32_t ***" );
   runTest< int32_t >( storage );

   WALBERLA_LOG_INFO_ON_ROOT( "*** Value Type = int64_t ***" );
   runTest< int64_t >( storage );

   return EXIT_SUCCESS;
}
