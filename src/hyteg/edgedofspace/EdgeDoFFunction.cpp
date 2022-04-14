/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "EdgeDoFFunction.hpp"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/edgedofspace/EdgeDoFAdditivePackInfo.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/edgedofspace/EdgeDoFPackInfo.hpp"
#include "hyteg/edgedofspace/generatedKernels/add_2D_macroface_edgedof_1_rhsfunction.hpp"
#include "hyteg/edgedofspace/generatedKernels/add_2D_macroface_edgedof_2_rhsfunctions.hpp"
#include "hyteg/edgedofspace/generatedKernels/add_2D_macroface_edgedof_3_rhsfunctions.hpp"
#include "hyteg/edgedofspace/generatedKernels/assign_2D_macroface_edgedof_1_rhsfunction.hpp"
#include "hyteg/edgedofspace/generatedKernels/assign_2D_macroface_edgedof_2_rhsfunctions.hpp"
#include "hyteg/edgedofspace/generatedKernels/assign_2D_macroface_edgedof_3_rhsfunctions.hpp"
#include "hyteg/edgedofspace/generatedKernels/assign_3D_macrocell_edgedof_1_rhsfunction.hpp"
#include "hyteg/edgedofspace/generatedKernels/assign_3D_macrocell_edgedof_2_rhsfunctions.hpp"
#include "hyteg/edgedofspace/generatedKernels/assign_3D_macrocell_edgedof_3_rhsfunctions.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/primitives/all.hpp"

#include "core/OpenMP.h"

namespace hyteg {

using walberla::int_c;

/// dummy function
template < typename ValueType >
EdgeDoFFunction< ValueType >::EdgeDoFFunction( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage )
: Function< EdgeDoFFunction< ValueType > >( name, storage )
, vertexDataID_( storage->generateInvalidPrimitiveDataID< MemoryDataHandling< FunctionMemory< ValueType >, Vertex >, Vertex >() )
, edgeDataID_( storage->generateInvalidPrimitiveDataID< MemoryDataHandling< FunctionMemory< ValueType >, Edge >, Edge >() )
, faceDataID_( storage->generateInvalidPrimitiveDataID< MemoryDataHandling< FunctionMemory< ValueType >, Face >, Face >() )
{}

template < typename ValueType >
EdgeDoFFunction< ValueType >::EdgeDoFFunction( const std::string&                         name,
                                               const std::shared_ptr< PrimitiveStorage >& storage,
                                               const uint_t&                              minLevel,
                                               const uint_t&                              maxLevel )
: EdgeDoFFunction( name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
{}

template < typename ValueType >
EdgeDoFFunction< ValueType >::EdgeDoFFunction( const std::string&                         name,
                                               const std::shared_ptr< PrimitiveStorage >& storage,
                                               const uint_t&                              minLevel,
                                               const uint_t&                              maxLevel,
                                               const BoundaryCondition&                   boundaryCondition )
: Function< EdgeDoFFunction< ValueType > >( name, storage, minLevel, maxLevel )
, boundaryCondition_( boundaryCondition )
{
   std::shared_ptr< MemoryDataHandling< FunctionMemory< ValueType >, Vertex > > vertexDataHandling =
       std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Vertex > >();

   std::shared_ptr< MemoryDataHandling< FunctionMemory< ValueType >, Edge > > edgeDataHandling =
       std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Edge > >();

   std::shared_ptr< MemoryDataHandling< FunctionMemory< ValueType >, Face > > faceDataHandling =
       std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Face > >();

   std::shared_ptr< MemoryDataHandling< FunctionMemory< ValueType >, Cell > > cellDataHandling =
       std::make_shared< MemoryDataHandling< FunctionMemory< ValueType >, Cell > >();

   storage->addVertexData( vertexDataID_, vertexDataHandling, name );
   storage->addEdgeData( edgeDataID_, edgeDataHandling, name );
   storage->addFaceData( faceDataID_, faceDataHandling, name );
   storage->addCellData( cellDataID_, cellDataHandling, name );

   for ( uint_t level = minLevel; level <= maxLevel; ++level )
   {
      for ( const auto & it : storage->getVertices() )
      {
         allocateMemory( level, *it.second );
      }
      for ( const auto & it : storage->getEdges() )
      {
         allocateMemory( level, *it.second );
      }
      for ( const auto & it : storage->getFaces() )
      {
         allocateMemory( level, *it.second );
      }
      for ( const auto & it : storage->getCells() )
      {
         allocateMemory( level, *it.second );
      }

      communicators_[level]->addPackInfo( std::make_shared< EdgeDoFPackInfo< ValueType > >(
          level, vertexDataID_, edgeDataID_, faceDataID_, cellDataID_, this->getStorage() ) );
      additiveCommunicators_[level]->addPackInfo( std::make_shared< EdgeDoFAdditivePackInfo< ValueType > >(
          level, vertexDataID_, edgeDataID_, faceDataID_, cellDataID_, this->getStorage() ) );
   }
}

template < typename ValueType >
bool EdgeDoFFunction< ValueType >::hasMemoryAllocated( const uint_t & level, const Vertex & vertex ) const
{
   WALBERLA_CHECK( this->getStorage()->vertexExistsLocally( vertex.getID() ) );
   return vertex.hasData( getVertexDataID() ) && vertex.getData( getVertexDataID() )->hasLevel( level );
}

template < typename ValueType >
bool EdgeDoFFunction< ValueType >::hasMemoryAllocated( const uint_t & level, const Edge & edge ) const
{
   WALBERLA_CHECK( this->getStorage()->edgeExistsLocally( edge.getID() ) );
   return edge.hasData( getEdgeDataID() ) && edge.getData( getEdgeDataID() )->hasLevel( level );
}

template < typename ValueType >
bool EdgeDoFFunction< ValueType >::hasMemoryAllocated( const uint_t & level, const Face & face ) const
{
   WALBERLA_CHECK( this->getStorage()->faceExistsLocally( face.getID() ) );
   return face.hasData( getFaceDataID() ) && face.getData( getFaceDataID() )->hasLevel( level );
}

template < typename ValueType >
bool EdgeDoFFunction< ValueType >::hasMemoryAllocated( const uint_t & level, const Cell & cell ) const
{
   WALBERLA_CHECK( this->getStorage()->cellExistsLocally( cell.getID() ) );
   return cell.hasData( getCellDataID() ) && cell.getData( getCellDataID() )->hasLevel( level );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::allocateMemory( const uint_t & level, const Vertex & vertex )
{
   WALBERLA_CHECK( this->getStorage()->vertexExistsLocally( vertex.getID() ) );
   WALBERLA_CHECK( vertex.hasData( getVertexDataID() ) )
   if ( hasMemoryAllocated( level, vertex ) )
      return;
   vertex.getData( getVertexDataID() )->addData( level, edgedof::edgeDoFMacroVertexFunctionMemorySize( level, vertex ), 0 );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::allocateMemory( const uint_t & level, const Edge & edge )
{
   WALBERLA_CHECK( this->getStorage()->edgeExistsLocally( edge.getID() ) );
   WALBERLA_CHECK( edge.hasData( getEdgeDataID() ) )
   if ( hasMemoryAllocated( level, edge ) )
      return;
   edge.getData( getEdgeDataID() )->addData( level, edgedof::edgeDoFMacroEdgeFunctionMemorySize( level, edge ), 0 );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::allocateMemory( const uint_t & level, const Face & face )
{
   WALBERLA_CHECK( this->getStorage()->faceExistsLocally( face.getID() ) );
   WALBERLA_CHECK( face.hasData( getFaceDataID() ) )
   if ( hasMemoryAllocated( level, face ) )
      return;
   face.getData( getFaceDataID() )->addData( level, edgedof::edgeDoFMacroFaceFunctionMemorySize( level, face ), 0 );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::allocateMemory( const uint_t & level, const Cell & cell )
{
   WALBERLA_CHECK( this->getStorage()->cellExistsLocally( cell.getID() ) );
   WALBERLA_CHECK( cell.hasData( getCellDataID() ) )
   if ( hasMemoryAllocated( level, cell ) )
      return;
   cell.getData( getCellDataID() )->addData( level, edgedof::edgeDoFMacroCellFunctionMemorySize( level, cell ), 0 );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::deleteMemory( const uint_t & level, const Vertex & vertex )
{
   WALBERLA_CHECK( this->getStorage()->vertexExistsLocally( vertex.getID() ) );
   if ( !hasMemoryAllocated( level, vertex ) )
      return;
   vertex.getData( getVertexDataID() )->deleteData( level );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::deleteMemory( const uint_t & level, const Edge & edge )
{
   WALBERLA_CHECK( this->getStorage()->edgeExistsLocally( edge.getID() ) );
   if ( !hasMemoryAllocated( level, edge ) )
      return;
   edge.getData( getEdgeDataID() )->deleteData( level );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::deleteMemory( const uint_t & level, const Face & face )
{
   WALBERLA_CHECK( this->getStorage()->faceExistsLocally( face.getID() ) );
   if ( !hasMemoryAllocated( level, face ) )
      return;
   face.getData( getFaceDataID() )->deleteData( level );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::deleteMemory( const uint_t & level, const Cell & cell )
{
   WALBERLA_CHECK( this->getStorage()->cellExistsLocally( cell.getID() ) );
   if ( !hasMemoryAllocated( level, cell ) )
      return;
   cell.getData( getCellDataID() )->deleteData( level );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::interpolate( ValueType constant, uint_t level, DoFType flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );
   interpolateByPrimitiveType< Edge >( constant, level, flag );
   interpolateByPrimitiveType< Face >( constant, level, flag );
   interpolateByPrimitiveType< Cell >( constant, level, flag );
   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::interpolate( ValueType constant, uint_t level, BoundaryUID boundaryUID ) const
{
  interpolate( [constant]( const Point3D& ){ return constant; }, level, boundaryUID );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                                uint_t                                              level,
                                                DoFType                                             flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) > exprExtended =
       [&expr]( const hyteg::Point3D& x, const std::vector< ValueType >& ) { return expr( x ); };
   interpolateExtended( exprExtended, {}, level, flag );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::interpolate( const std::function< ValueType( const Point3D& ) >& expr,
                                                uint_t                                              level,
                                                BoundaryUID                                         boundaryUID ) const
{
   if ( isDummy() )
   {
      return;
   }
   std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) > exprExtended =
       [&expr]( const hyteg::Point3D& x, const std::vector< ValueType >& ) { return expr( x ); };
   interpolateExtended( exprExtended, {}, level, boundaryUID );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::interpolateExtended(
    const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
    const std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >&   srcFunctions,
    uint_t                                                                               level,
    DoFType                                                                              flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > > srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > srcCellIDs;

   for ( const EdgeDoFFunction& function : srcFunctions )
   {
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         edgedof::macroedge::interpolate< ValueType >( level, edge, edgeDataID_, srcEdgeIDs, expr );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c(i)] );   

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         edgedof::macroface::interpolate< ValueType >( level, face, faceDataID_, srcFaceIDs, expr );
      }
   }

   if ( level >= 1 )
   {
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c(i)] );

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            edgedof::macrocell::interpolate< ValueType >( level, cell, cellDataID_, srcCellIDs, expr );
         }
      }
   }
   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::interpolateExtended(
    const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
    const std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >&   srcFunctions,
    uint_t                                                                               level,
    BoundaryUID                                                                          boundaryUID ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > > srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > srcCellIDs;

   for ( const EdgeDoFFunction& function : srcFunctions )
   {
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c(i)] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( edge.getMeshBoundaryFlag() ) == boundaryUID )
      {
         edgedof::macroedge::interpolate< ValueType >( level, edge, edgeDataID_, srcEdgeIDs, expr );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c(i)] );

      if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( face.getMeshBoundaryFlag() ) == boundaryUID )
      {
         edgedof::macroface::interpolate< ValueType >( level, face, faceDataID_, srcFaceIDs, expr );
      }
   }

   if ( level >= 1 )
   {
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c(i)] );

         if ( boundaryCondition_.getBoundaryUIDFromMeshFlag( cell.getMeshBoundaryFlag() ) == boundaryUID )
         {
            edgedof::macrocell::interpolate< ValueType >( level, cell, cellDataID_, srcCellIDs, expr );
         }
      }
   }
   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::setToZero( uint_t level ) const
{
   if ( isDummy() )
   {
      return;
   }

   this->startTiming( "setToZero" );

   for ( const auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;
      vertex.getData( vertexDataID_ )->setToZero( level );
   }

   for ( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;
      edge.getData( edgeDataID_ )->setToZero( level );
   }

   for ( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;
      face.getData( faceDataID_ )->setToZero( level );
   }

   for ( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      cell.getData( cellDataID_ )->setToZero( level );
   }

   this->stopTiming( "setToZero" );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::swap( const EdgeDoFFunction< ValueType >& other,
                                         const uint_t&                       level,
                                         const DoFType&                      flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Swap" );

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         edgedof::macroedge::swap< ValueType >( level, edge, other.getEdgeDataID(), edgeDataID_ );
      }
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         edgedof::macroface::swap< ValueType >( level, face, other.getFaceDataID(), faceDataID_ );
      }
   }

   if ( level >= 1 )
   {
      for ( auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            edgedof::macrocell::swap< ValueType >( level, cell, other.getCellDataID(), cellDataID_ );
         }
      }
   }

   this->stopTiming( "Swap" );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::copyFrom( const EdgeDoFFunction< ValueType >& other, const uint_t& level ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Copy" );

   for ( auto& it : this->getStorage()->getVertices() )
   {
      auto primitiveID = it.first;
      WALBERLA_ASSERT( other.getStorage()->vertexExistsLocally( primitiveID ) )
      this->getStorage()
          ->getVertex( primitiveID )
          ->getData( vertexDataID_ )
          ->copyFrom( *other.getStorage()->getVertex( primitiveID )->getData( other.getVertexDataID() ), level );
   }

   for ( auto& it : this->getStorage()->getEdges() )
   {
      auto primitiveID = it.first;
      WALBERLA_ASSERT( other.getStorage()->edgeExistsLocally( primitiveID ) )
      this->getStorage()
          ->getEdge( primitiveID )
          ->getData( edgeDataID_ )
          ->copyFrom( *other.getStorage()->getEdge( primitiveID )->getData( other.getEdgeDataID() ), level );
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      auto primitiveID = it.first;
      WALBERLA_ASSERT( other.getStorage()->faceExistsLocally( primitiveID ) )
      this->getStorage()
          ->getFace( primitiveID )
          ->getData( faceDataID_ )
          ->copyFrom( *other.getStorage()->getFace( primitiveID )->getData( other.getFaceDataID() ), level );
   }

   for ( auto& it : this->getStorage()->getCells() )
   {
      auto primitiveID = it.first;
      WALBERLA_ASSERT( other.getStorage()->cellExistsLocally( primitiveID ) )
      this->getStorage()
          ->getCell( primitiveID )
          ->getData( cellDataID_ )
          ->copyFrom( *other.getStorage()->getCell( primitiveID )->getData( other.getCellDataID() ), level );
   }

   this->stopTiming( "Copy" );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::copyFrom( const EdgeDoFFunction< ValueType >&            other,
                                             const uint_t&                                  level,
                                             const std::map< PrimitiveID::IDType, uint_t >& localPrimitiveIDsToRank,
                                             const std::map< PrimitiveID::IDType, uint_t >& otherPrimitiveIDsToRank ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Copy" );

   walberla::mpi::BufferSystem bufferSystem( walberla::mpi::MPIManager::instance()->comm(), 9563 );
   std::set< walberla::mpi::MPIRank > receiverRanks;
   for ( auto it : localPrimitiveIDsToRank )
   {
      receiverRanks.insert( walberla::mpi::MPIRank( it.second ) );
   }
   bufferSystem.setReceiverInfo( receiverRanks, true );

   for ( auto& it : other.getStorage()->getVertices() )
   {
      PrimitiveID::IDType otherPrimitiveID = it.first;
      WALBERLA_CHECK_GREATER( otherPrimitiveIDsToRank.count( otherPrimitiveID ), 0 );
      auto otherData     = it.second->getData( other.getVertexDataID() )->getPointer( level );
      auto otherDataSize = it.second->getData( other.getVertexDataID() )->getSize( level );
      auto targetRank    = otherPrimitiveIDsToRank.at( otherPrimitiveID );
      bufferSystem.sendBuffer( targetRank ) << otherPrimitiveID;
      bufferSystem.sendBuffer( targetRank ) << uint_c( 0 );
      bufferSystem.sendBuffer( targetRank ) << otherDataSize;
      for ( uint_t i = 0; i < otherDataSize; i++ )
         bufferSystem.sendBuffer( targetRank ) << otherData[i];
   }

   for ( auto& it : other.getStorage()->getEdges() )
   {
      PrimitiveID::IDType otherPrimitiveID = it.first;
      WALBERLA_CHECK_GREATER( otherPrimitiveIDsToRank.count( otherPrimitiveID ), 0 );
      auto otherData     = it.second->getData( other.getEdgeDataID() )->getPointer( level );
      auto otherDataSize = it.second->getData( other.getEdgeDataID() )->getSize( level );
      auto targetRank    = otherPrimitiveIDsToRank.at( otherPrimitiveID );
      bufferSystem.sendBuffer( targetRank ) << otherPrimitiveID;
      bufferSystem.sendBuffer( targetRank ) << uint_c( 1 );
      bufferSystem.sendBuffer( targetRank ) << otherDataSize;
      for ( uint_t i = 0; i < otherDataSize; i++ )
         bufferSystem.sendBuffer( targetRank ) << otherData[i];
   }

   for ( auto& it : other.getStorage()->getFaces() )
   {
      PrimitiveID::IDType otherPrimitiveID = it.first;
      WALBERLA_CHECK_GREATER( otherPrimitiveIDsToRank.count( otherPrimitiveID ), 0 );
      auto otherData     = it.second->getData( other.getFaceDataID() )->getPointer( level );
      auto otherDataSize = it.second->getData( other.getFaceDataID() )->getSize( level );
      auto targetRank    = otherPrimitiveIDsToRank.at( otherPrimitiveID );
      bufferSystem.sendBuffer( targetRank ) << otherPrimitiveID;
      bufferSystem.sendBuffer( targetRank ) << uint_c( 2 );
      bufferSystem.sendBuffer( targetRank ) << otherDataSize;
      for ( uint_t i = 0; i < otherDataSize; i++ )
         bufferSystem.sendBuffer( targetRank ) << otherData[i];
   }

   for ( auto& it : other.getStorage()->getCells() )
   {
      PrimitiveID::IDType otherPrimitiveID = it.first;
      WALBERLA_CHECK_GREATER( otherPrimitiveIDsToRank.count( otherPrimitiveID ), 0 );
      auto otherData     = it.second->getData( other.getCellDataID() )->getPointer( level );
      auto otherDataSize = it.second->getData( other.getCellDataID() )->getSize( level );
      auto targetRank    = otherPrimitiveIDsToRank.at( otherPrimitiveID );
      bufferSystem.sendBuffer( targetRank ) << otherPrimitiveID;
      bufferSystem.sendBuffer( targetRank ) << uint_c( 3 );
      bufferSystem.sendBuffer( targetRank ) << otherDataSize;
      for ( uint_t i = 0; i < otherDataSize; i++ )
         bufferSystem.sendBuffer( targetRank ) << otherData[i];
   }

   bufferSystem.sendAll();

   for ( auto pkg = bufferSystem.begin(); pkg != bufferSystem.end(); ++pkg )
   {
      while ( !pkg.buffer().isEmpty() )
      {
         PrimitiveID::IDType otherID;
         uint_t              primitiveType = 4;
         uint_t              dataSize      = 0;
         ValueType           value;
         ValueType*          dstPointer;

         pkg.buffer() >> otherID;
         pkg.buffer() >> primitiveType;
         pkg.buffer() >> dataSize;

         WALBERLA_CHECK( this->getStorage()->primitiveExistsLocally( PrimitiveID( otherID ) ) );

         switch ( primitiveType )
         {
         case 0:
            dstPointer = this->getStorage()->getVertex( PrimitiveID( otherID ) )->getData( vertexDataID_ )->getPointer( level );
            break;
         case 1:
            dstPointer = this->getStorage()->getEdge( PrimitiveID( otherID ) )->getData( edgeDataID_ )->getPointer( level );
            break;
         case 2:
            dstPointer = this->getStorage()->getFace( PrimitiveID( otherID ) )->getData( faceDataID_ )->getPointer( level );
            break;
         case 3:
            dstPointer = this->getStorage()->getCell( PrimitiveID( otherID ) )->getData( cellDataID_ )->getPointer( level );
            break;
         default:
         WALBERLA_ABORT( "Invalid primitive type" )
         }

         for ( uint_t i = 0; i < dataSize; i++ )
         {
            pkg.buffer() >> value;
            dstPointer[i] = value;
         }
      }
   }

   this->stopTiming( "Copy" );
}


template < typename ValueType >
void macroFaceAssign( const uint_t&                                                              level,
                      Face&                                                                      face,
                      const std::vector< ValueType >&                                            scalars,
                      const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcFaceIDs,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstFaceID )
{
   edgedof::macroface::assign< ValueType >( level, face, scalars, srcFaceIDs, dstFaceID );
}

template <>
void macroFaceAssign< double >( const uint_t&                                                           level,
                                Face&                                                                   face,
                                const std::vector< double >&                                            scalars,
                                const std::vector< PrimitiveDataID< FunctionMemory< double >, Face > >& srcFaceIDs,
                                const PrimitiveDataID< FunctionMemory< double >, Face >&                dstFaceID )
{
   typedef edgedof::EdgeDoFOrientation eo;
   std::map< eo, uint_t >              firstIdx;
   for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
      firstIdx[e] = edgedof::macroface::index( level, 0, 0, e );

   if ( globalDefines::useGeneratedKernels && scalars.size() == 1 )
   {
      auto dstData = face.getData( dstFaceID )->getPointer( level );
      auto srcData = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto scalar  = scalars.at( 0 );
      edgedof::macroface::generated::assign_2D_macroface_edgedof_1_rhs_function( &dstData[firstIdx[eo::X]],
                                                                                 &dstData[firstIdx[eo::XY]],
                                                                                 &dstData[firstIdx[eo::Y]],
                                                                                 &srcData[firstIdx[eo::X]],
                                                                                 &srcData[firstIdx[eo::XY]],
                                                                                 &srcData[firstIdx[eo::Y]],
                                                                                 scalar,
                                                                                 static_cast< int32_t >( level ) );
   }
   else if ( globalDefines::useGeneratedKernels && scalars.size() == 2 )
   {
      auto dstData  = face.getData( dstFaceID )->getPointer( level );
      auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      edgedof::macroface::generated::assign_2D_macroface_edgedof_2_rhs_functions( &dstData[firstIdx[eo::X]],
                                                                                  &dstData[firstIdx[eo::XY]],
                                                                                  &dstData[firstIdx[eo::Y]],
                                                                                  &srcData0[firstIdx[eo::X]],
                                                                                  &srcData0[firstIdx[eo::XY]],
                                                                                  &srcData0[firstIdx[eo::Y]],
                                                                                  &srcData1[firstIdx[eo::X]],
                                                                                  &srcData1[firstIdx[eo::XY]],
                                                                                  &srcData1[firstIdx[eo::Y]],
                                                                                  scalar0,
                                                                                  scalar1,
                                                                                  static_cast< int32_t >( level ) );
   }
   else if ( globalDefines::useGeneratedKernels && scalars.size() == 3 )
   {
      auto dstData  = face.getData( dstFaceID )->getPointer( level );
      auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
      auto srcData2 = face.getData( srcFaceIDs.at( 2 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      auto scalar2  = scalars.at( 2 );
      edgedof::macroface::generated::assign_2D_macroface_edgedof_3_rhs_functions( &dstData[firstIdx[eo::X]],
                                                                                  &dstData[firstIdx[eo::XY]],
                                                                                  &dstData[firstIdx[eo::Y]],
                                                                                  &srcData0[firstIdx[eo::X]],
                                                                                  &srcData0[firstIdx[eo::XY]],
                                                                                  &srcData0[firstIdx[eo::Y]],
                                                                                  &srcData1[firstIdx[eo::X]],
                                                                                  &srcData1[firstIdx[eo::XY]],
                                                                                  &srcData1[firstIdx[eo::Y]],
                                                                                  &srcData2[firstIdx[eo::X]],
                                                                                  &srcData2[firstIdx[eo::XY]],
                                                                                  &srcData2[firstIdx[eo::Y]],
                                                                                  scalar0,
                                                                                  scalar1,
                                                                                  scalar2,
                                                                                  static_cast< int32_t >( level ) );
   }
   else
   {
      edgedof::macroface::assign< double >( level, face, scalars, srcFaceIDs, dstFaceID );
   }
}

template < typename ValueType >
void macroCellAssign( const uint_t&                                                              level,
                      Cell&                                                                      cell,
                      const std::vector< ValueType >&                                            scalars,
                      const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > >& srcCellIDs,
                      const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&                dstCellID )
{
   edgedof::macrocell::assign< ValueType >( level, cell, scalars, srcCellIDs, dstCellID );
}

template <>
void macroCellAssign< double >( const uint_t&                                                           level,
                                Cell&                                                                   cell,
                                const std::vector< double >&                                            scalars,
                                const std::vector< PrimitiveDataID< FunctionMemory< double >, Cell > >& srcCellIDs,
                                const PrimitiveDataID< FunctionMemory< double >, Cell >&                dstCellID )
{
   typedef edgedof::EdgeDoFOrientation eo;
   if ( globalDefines::useGeneratedKernels && scalars.size() == 1 )
   {
      auto                   dstData = cell.getData( dstCellID )->getPointer( level );
      auto                   srcData = cell.getData( srcCellIDs.at( 0 ) )->getPointer( level );
      auto                   scalar  = scalars.at( 0 );
      std::map< eo, uint_t > firstIdx;
      for ( auto e : edgedof::allEdgeDoFOrientations )
         firstIdx[e] = edgedof::macrocell::index( level, 0, 0, 0, e );
      edgedof::macrocell::generated::assign_3D_macrocell_edgedof_1_rhs_function( &dstData[firstIdx[eo::X]],
                                                                                 &dstData[firstIdx[eo::XY]],
                                                                                 &dstData[firstIdx[eo::XYZ]],
                                                                                 &dstData[firstIdx[eo::XZ]],
                                                                                 &dstData[firstIdx[eo::Y]],
                                                                                 &dstData[firstIdx[eo::YZ]],
                                                                                 &dstData[firstIdx[eo::Z]],
                                                                                 &srcData[firstIdx[eo::X]],
                                                                                 &srcData[firstIdx[eo::XY]],
                                                                                 &srcData[firstIdx[eo::XYZ]],
                                                                                 &srcData[firstIdx[eo::XZ]],
                                                                                 &srcData[firstIdx[eo::Y]],
                                                                                 &srcData[firstIdx[eo::YZ]],
                                                                                 &srcData[firstIdx[eo::Z]],
                                                                                 scalar,
                                                                                 static_cast< int32_t >( level ) );
   }
   else if ( globalDefines::useGeneratedKernels && scalars.size() == 2 )
   {
      auto                   dstData  = cell.getData( dstCellID )->getPointer( level );
      auto                   srcData0 = cell.getData( srcCellIDs.at( 0 ) )->getPointer( level );
      auto                   srcData1 = cell.getData( srcCellIDs.at( 1 ) )->getPointer( level );
      auto                   scalar0  = scalars.at( 0 );
      auto                   scalar1  = scalars.at( 1 );
      std::map< eo, uint_t > firstIdx;
      for ( auto e : edgedof::allEdgeDoFOrientations )
         firstIdx[e] = edgedof::macrocell::index( level, 0, 0, 0, e );

      edgedof::macrocell::generated::assign_3D_macrocell_edgedof_2_rhs_functions( &dstData[firstIdx[eo::X]],
                                                                                  &dstData[firstIdx[eo::XY]],
                                                                                  &dstData[firstIdx[eo::XYZ]],
                                                                                  &dstData[firstIdx[eo::XZ]],
                                                                                  &dstData[firstIdx[eo::Y]],
                                                                                  &dstData[firstIdx[eo::YZ]],
                                                                                  &dstData[firstIdx[eo::Z]],
                                                                                  &srcData0[firstIdx[eo::X]],
                                                                                  &srcData0[firstIdx[eo::XY]],
                                                                                  &srcData0[firstIdx[eo::XYZ]],
                                                                                  &srcData0[firstIdx[eo::XZ]],
                                                                                  &srcData0[firstIdx[eo::Y]],
                                                                                  &srcData0[firstIdx[eo::YZ]],
                                                                                  &srcData0[firstIdx[eo::Z]],
                                                                                  &srcData1[firstIdx[eo::X]],
                                                                                  &srcData1[firstIdx[eo::XY]],
                                                                                  &srcData1[firstIdx[eo::XYZ]],
                                                                                  &srcData1[firstIdx[eo::XZ]],
                                                                                  &srcData1[firstIdx[eo::Y]],
                                                                                  &srcData1[firstIdx[eo::YZ]],
                                                                                  &srcData1[firstIdx[eo::Z]],
                                                                                  scalar0,
                                                                                  scalar1,
                                                                                  static_cast< int32_t >( level ) );
   }
   else if ( globalDefines::useGeneratedKernels && scalars.size() == 3 )
   {
      auto                   dstData  = cell.getData( dstCellID )->getPointer( level );
      auto                   srcData0 = cell.getData( srcCellIDs.at( 0 ) )->getPointer( level );
      auto                   srcData1 = cell.getData( srcCellIDs.at( 1 ) )->getPointer( level );
      auto                   srcData2 = cell.getData( srcCellIDs.at( 2 ) )->getPointer( level );
      auto                   scalar0  = scalars.at( 0 );
      auto                   scalar1  = scalars.at( 1 );
      auto                   scalar2  = scalars.at( 2 );
      std::map< eo, uint_t > firstIdx;
      for ( auto e : edgedof::allEdgeDoFOrientations )
         firstIdx[e] = edgedof::macrocell::index( level, 0, 0, 0, e );

      edgedof::macrocell::generated::assign_3D_macrocell_edgedof_3_rhs_functions( &dstData[firstIdx[eo::X]],
                                                                                  &dstData[firstIdx[eo::XY]],
                                                                                  &dstData[firstIdx[eo::XYZ]],
                                                                                  &dstData[firstIdx[eo::XZ]],
                                                                                  &dstData[firstIdx[eo::Y]],
                                                                                  &dstData[firstIdx[eo::YZ]],
                                                                                  &dstData[firstIdx[eo::Z]],
                                                                                  &srcData0[firstIdx[eo::X]],
                                                                                  &srcData0[firstIdx[eo::XY]],
                                                                                  &srcData0[firstIdx[eo::XYZ]],
                                                                                  &srcData0[firstIdx[eo::XZ]],
                                                                                  &srcData0[firstIdx[eo::Y]],
                                                                                  &srcData0[firstIdx[eo::YZ]],
                                                                                  &srcData0[firstIdx[eo::Z]],
                                                                                  &srcData1[firstIdx[eo::X]],
                                                                                  &srcData1[firstIdx[eo::XY]],
                                                                                  &srcData1[firstIdx[eo::XYZ]],
                                                                                  &srcData1[firstIdx[eo::XZ]],
                                                                                  &srcData1[firstIdx[eo::Y]],
                                                                                  &srcData1[firstIdx[eo::YZ]],
                                                                                  &srcData1[firstIdx[eo::Z]],
                                                                                  &srcData2[firstIdx[eo::X]],
                                                                                  &srcData2[firstIdx[eo::XY]],
                                                                                  &srcData2[firstIdx[eo::XYZ]],
                                                                                  &srcData2[firstIdx[eo::XZ]],
                                                                                  &srcData2[firstIdx[eo::Y]],
                                                                                  &srcData2[firstIdx[eo::YZ]],
                                                                                  &srcData2[firstIdx[eo::Z]],
                                                                                  scalar0,
                                                                                  scalar1,
                                                                                  scalar2,
                                                                                  static_cast< int32_t >( level ) );
   }
   else
   {
      edgedof::macrocell::assign< double >( level, cell, scalars, srcCellIDs, dstCellID );
   }
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::assign(
    const std::vector< ValueType >&                                                    scalars,
    const std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >& functions,
    size_t                                                                             level,
    DoFType                                                                            flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Assign" );
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > > srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > srcCellIDs;

   for ( const EdgeDoFFunction< ValueType >& function : functions )
   {
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         edgedof::macroedge::assign< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         macroFaceAssign< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_ );
      }
   }

   if ( level >= 1 )
   {
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
      #pragma omp parallel for default(shared)
      #endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c(i)] );

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            macroCellAssign< ValueType >( level, cell, scalars, srcCellIDs, cellDataID_ );
         }
      }
   }

   this->stopTiming( "Assign" );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::add( ValueType scalar, uint_t level, DoFType flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Add (scalar)" );

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         edgedof::macroedge::add< ValueType >( level, edge, scalar, edgeDataID_ );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         edgedof::macroface::add< ValueType >( level, face, scalar, faceDataID_ );
      }
   }

   if ( level >= 1 )
   {
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c(i)] );

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            edgedof::macrocell::add< ValueType >( level, cell, scalar, cellDataID_ );
         }
      }
   }

   this->stopTiming( "Add (scalar)" );
}

template < typename ValueType >
void macroFaceAdd( const uint_t&                                                              level,
                   Face&                                                                      face,
                   const std::vector< ValueType >&                                            scalars,
                   const std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >& srcFaceIDs,
                   const PrimitiveDataID< FunctionMemory< ValueType >, Face >&                dstFaceID )
{
   edgedof::macroface::add< ValueType >( level, face, scalars, srcFaceIDs, dstFaceID );
}

template <>
void macroFaceAdd< double >( const uint_t&                                                           level,
                             Face&                                                                   face,
                             const std::vector< double >&                                            scalars,
                             const std::vector< PrimitiveDataID< FunctionMemory< double >, Face > >& srcFaceIDs,
                             const PrimitiveDataID< FunctionMemory< double >, Face >&                dstFaceID )
{
   typedef edgedof::EdgeDoFOrientation eo;
   std::map< eo, uint_t >              firstIdx;
   for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
      firstIdx[e] = edgedof::macroface::index( level, 0, 0, e );

   if ( globalDefines::useGeneratedKernels && scalars.size() == 1 )
   {
      auto dstData = face.getData( dstFaceID )->getPointer( level );
      auto srcData = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto scalar  = scalars.at( 0 );
      edgedof::macroface::generated::add_2D_macroface_edgedof_1_rhs_function( &dstData[firstIdx[eo::X]],
                                                                              &dstData[firstIdx[eo::XY]],
                                                                              &dstData[firstIdx[eo::Y]],
                                                                              &srcData[firstIdx[eo::X]],
                                                                              &srcData[firstIdx[eo::XY]],
                                                                              &srcData[firstIdx[eo::Y]],
                                                                              scalar,
                                                                              static_cast< int32_t >( level ) );
   }
   else if ( globalDefines::useGeneratedKernels && scalars.size() == 2 )
   {
      auto dstData  = face.getData( dstFaceID )->getPointer( level );
      auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      edgedof::macroface::generated::add_2D_macroface_edgedof_2_rhs_functions( &dstData[firstIdx[eo::X]],
                                                                               &dstData[firstIdx[eo::XY]],
                                                                               &dstData[firstIdx[eo::Y]],
                                                                               &srcData0[firstIdx[eo::X]],
                                                                               &srcData0[firstIdx[eo::XY]],
                                                                               &srcData0[firstIdx[eo::Y]],
                                                                               &srcData1[firstIdx[eo::X]],
                                                                               &srcData1[firstIdx[eo::XY]],
                                                                               &srcData1[firstIdx[eo::Y]],
                                                                               scalar0,
                                                                               scalar1,
                                                                               static_cast< int32_t >( level ) );
   }
   else if ( globalDefines::useGeneratedKernels && scalars.size() == 3 )
   {
      auto dstData  = face.getData( dstFaceID )->getPointer( level );
      auto srcData0 = face.getData( srcFaceIDs.at( 0 ) )->getPointer( level );
      auto srcData1 = face.getData( srcFaceIDs.at( 1 ) )->getPointer( level );
      auto srcData2 = face.getData( srcFaceIDs.at( 2 ) )->getPointer( level );
      auto scalar0  = scalars.at( 0 );
      auto scalar1  = scalars.at( 1 );
      auto scalar2  = scalars.at( 2 );
      edgedof::macroface::generated::add_2D_macroface_edgedof_3_rhs_functions( &dstData[firstIdx[eo::X]],
                                                                               &dstData[firstIdx[eo::XY]],
                                                                               &dstData[firstIdx[eo::Y]],
                                                                               &srcData0[firstIdx[eo::X]],
                                                                               &srcData0[firstIdx[eo::XY]],
                                                                               &srcData0[firstIdx[eo::Y]],
                                                                               &srcData1[firstIdx[eo::X]],
                                                                               &srcData1[firstIdx[eo::XY]],
                                                                               &srcData1[firstIdx[eo::Y]],
                                                                               &srcData2[firstIdx[eo::X]],
                                                                               &srcData2[firstIdx[eo::XY]],
                                                                               &srcData2[firstIdx[eo::Y]],
                                                                               scalar0,
                                                                               scalar1,
                                                                               scalar2,
                                                                               static_cast< int32_t >( level ) );
   }
   else
   {
      edgedof::macroface::add< double >( level, face, scalars, srcFaceIDs, dstFaceID );
   }
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::add(
    const std::vector< ValueType >&                                                    scalars,
    const std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >& functions,
    size_t                                                                             level,
    DoFType                                                                            flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Add" );
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > > srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > srcCellIDs;

   for ( const EdgeDoFFunction< ValueType >& function : functions )
   {
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         edgedof::macroedge::add< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         macroFaceAdd< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_ );
      }
   }

   if ( level >= 1 )
   {
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
      #pragma omp parallel for default(shared)
      #endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c(i)] );

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            edgedof::macrocell::add< ValueType >( level, cell, scalars, srcCellIDs, cellDataID_ );
         }
      }
   }

   this->stopTiming( "Add" );
}

template < typename ValueType >
ValueType EdgeDoFFunction< ValueType >::dotLocal( const EdgeDoFFunction< ValueType >& rhs,
                                                  const uint_t                        level,
                                                  const DoFType                       flag ) const
{
   if ( isDummy() )
   {
      return ValueType( 0 );
   }
   this->startTiming( "Dot (local)" );
   auto scalarProduct = ValueType( 0 );

   ValueType scalarProductEdges = 0;
   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for reduction(+: scalarProductEdges)
   #endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[ uint_c(i) ] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         scalarProductEdges += edgedof::macroedge::dot< ValueType >( level, edge, edgeDataID_, rhs.edgeDataID_ );
      }
   }
   scalarProduct += scalarProductEdges;

   ValueType scalarProductFaces = 0;
   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for reduction(+: scalarProductFaces)
   #endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[ uint_c(i) ] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         scalarProductFaces += edgedof::macroface::dot< ValueType >( level, face, faceDataID_, rhs.faceDataID_ );
      }
   }
   scalarProduct += scalarProductFaces;

   ValueType scalarProductCells = 0;
   if ( level >= 1 )
   {
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
      #pragma omp parallel for reduction(+: scalarProductCells)
      #endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[ uint_c(i) ] );

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            scalarProductCells += edgedof::macrocell::dot< ValueType >( level, cell, cellDataID_, rhs.cellDataID_ );
         }
      }
   }
   scalarProduct += scalarProductCells;

   this->stopTiming( "Dot (local)" );

   return scalarProduct;
}

template < typename ValueType >
ValueType EdgeDoFFunction< ValueType >::sumGlobal( const uint_t& level, const DoFType& flag, const bool& absolute ) const
{
   ValueType sum = sumLocal( level, flag, absolute );
   this->startTiming( "Sum (reduce)" );
   walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
   this->stopTiming( "Sum (reduce)" );
   return sum;
}

template < typename ValueType >
ValueType EdgeDoFFunction< ValueType >::sumLocal( const uint_t& level, const DoFType& flag, const bool& absolute ) const
{
   if ( isDummy() )
   {
      return ValueType( 0 );
   }
   this->startTiming( "Sum (local)" );
   auto sum = ValueType( 0 );

   for ( const auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         sum += edgedof::macroedge::sum< ValueType >( level, edge, edgeDataID_, absolute );
      }
   }

   for ( const auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         sum += edgedof::macroface::sum< ValueType >( level, face, faceDataID_, absolute );
      }
   }

   for ( const auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         sum += edgedof::macrocell::sum< ValueType >( level, cell, cellDataID_, absolute );
      }
   }
   this->stopTiming( "Sum (local)" );
   return sum;
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::enumerate( uint_t level ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Enumerate" );

   auto counter = static_cast< ValueType >( hyteg::numberOfLocalDoFs< EdgeDoFFunctionTag >( *( this->getStorage() ), level ) );

   std::vector< ValueType > dofs_per_rank = walberla::mpi::allGather( counter );

   auto startOnRank = ValueType( 0 );

   for ( uint_t i = 0; i < uint_c( walberla::MPIManager::instance()->rank() ); ++i )
   {
      startOnRank += dofs_per_rank[i];
   }

   enumerate( level, startOnRank );
   this->stopTiming( "Enumerate" );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::enumerate( uint_t level, ValueType& offset ) const
{
   if ( isDummy() )
   {
      return;
   }

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;
      edgedof::macroedge::enumerate< ValueType >( level, edge, edgeDataID_, offset );
   }

   if ( level >= 1 )
   {
      for ( auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;
         edgedof::macroface::enumerate< ValueType >( level, face, faceDataID_, offset );
      }
   }

   if ( level >= 1 )
   {
      for ( auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;
         edgedof::macrocell::enumerate< ValueType >( level, cell, cellDataID_, offset );
      }
   }

   communication::syncFunctionBetweenPrimitives( *this, level );
}

template < typename ValueType >
ValueType EdgeDoFFunction< ValueType >::getMaxValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   if ( isDummy() )
   {
      return ValueType( 0 );
   }

   auto localMax = -std::numeric_limits< ValueType >::max();

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         localMax = std::max( localMax, edgedof::macroedge::getMaxValue< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         localMax = std::max( localMax, edgedof::macroface::getMaxValue< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      localMax   = std::max( localMax, edgedof::macrocell::getMaxValue< ValueType >( level, cell, cellDataID_ ) );
   }

   if ( mpiReduce )
   {
      walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );
   }

   return localMax;
}

template < typename ValueType >
ValueType EdgeDoFFunction< ValueType >::getMinValue( uint_t level, DoFType flag, bool mpiReduce ) const
{
   if ( isDummy() )
   {
      return ValueType( 0 );
   }

   auto localMin = std::numeric_limits< ValueType >::max();

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         localMin = std::min( localMin, edgedof::macroedge::getMinValue< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         localMin = std::min( localMin, edgedof::macroface::getMinValue< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      localMin   = std::min( localMin, edgedof::macrocell::getMinValue< ValueType >( level, cell, cellDataID_ ) );
   }

   if ( mpiReduce )
   {
      walberla::mpi::allReduceInplace( localMin, walberla::mpi::MIN, walberla::mpi::MPIManager::instance()->comm() );
   }

   return localMin;
}

template < typename ValueType >
ValueType EdgeDoFFunction< ValueType >::getMaxMagnitude( uint_t level, DoFType flag, bool mpiReduce ) const
{
   if ( isDummy() )
   {
      return ValueType( 0 );
   }
   auto localMax = ValueType( 0.0 );

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge&         edge   = *it.second;
      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         localMax = std::max( localMax, edgedof::macroedge::getMaxMagnitude< ValueType >( level, edge, edgeDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         localMax = std::max( localMax, edgedof::macroface::getMaxMagnitude< ValueType >( level, face, faceDataID_ ) );
      }
   }

   for ( auto& it : this->getStorage()->getCells() )
   {
      Cell& cell = *it.second;
      localMax   = std::max( localMax, edgedof::macrocell::getMaxMagnitude< ValueType >( level, cell, cellDataID_ ) );
   }

   if ( mpiReduce )
   {
      walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );
   }

   return localMax;
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::multElementwise(
    const std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >& functions,
    const uint_t                                                                       level,
    const DoFType                                                                      flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Multiply elementwise" );

   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > > srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > srcFaceIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Cell > > srcCellIDs;

   for ( const EdgeDoFFunction& function : functions )
   {
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
      srcCellIDs.push_back( function.cellDataID_ );
   }

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         edgedof::macroedge::multElementwise< ValueType >( level, edge, srcEdgeIDs, edgeDataID_ );
      }
   }

   std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *this->getStorage()->getFace( faceIDs[uint_c(i)] );

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         edgedof::macroface::multElementwise< ValueType >( level, face, srcFaceIDs, faceDataID_ );
      }
   }

   std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
   #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c(i)] );
      if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
      {
         edgedof::macrocell::multElementwise< ValueType >( level, cell, srcCellIDs, cellDataID_ );
      }
   }
   this->stopTiming( "Multiply elementwise" );
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::setLocalCommunicationMode(
    const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode )
{
   if ( isDummy() )
   {
      return;
   }
   for ( auto& communicator : communicators_ )
   {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
   }
   for ( auto& communicator : additiveCommunicators_ )
   {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
   }
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::invertElementwise( const uint_t level, const DoFType flag, bool workOnHalos ) const
{
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( flag );
   WALBERLA_UNUSED( workOnHalos );
   WALBERLA_ABORT( "EdgeDoFFunction< ValueType >::invertElementwise not available for requested ValueType" );
}

template < typename ValueType >
template < typename PrimitiveType >
void EdgeDoFFunction< ValueType >::interpolateByPrimitiveType( const ValueType& constant, uint_t level, DoFType flag ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Interpolate" );

   if ( std::is_same< PrimitiveType, Edge >::value )
   {
      std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
      for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
      {
         Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c(i)] );

         if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
         {
            edgedof::macroedge::interpolate( level, edge, edgeDataID_, constant );
         }
      }
   }
   else if ( std::is_same< PrimitiveType, Face >::value )
   {
      std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
      for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
      {
         Face& face = *this->getStorage()->getFace( faceIDs[uint_c(i)] );

         if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
         {
            edgedof::macroface::interpolate( level, face, faceDataID_, constant );
         }
      }
   }
   else if ( std::is_same< PrimitiveType, Cell >::value )
   {
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
      #ifdef WALBERLA_BUILD_WITH_OPENMP
   #pragma omp parallel for default(shared)
   #endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c(i)] );

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            edgedof::macrocell::interpolate( level, cell, cellDataID_, constant );
         }
      }
   }

   this->stopTiming( "Interpolate" );
}

uint_t edgedof::edgeDoFMacroVertexFunctionMemorySize( const uint_t& level, const Primitive& primitive )
{
   WALBERLA_UNUSED( level );
   return primitive.getNumNeighborEdges() + primitive.getNumNeighborFaces();
}

uint_t edgedof::edgeDoFMacroEdgeFunctionMemorySize( const uint_t& level, const Primitive& primitive )
{
   if ( level == 0 )
   {
      return 1 /* on edge */ + 2 * primitive.getNumNeighborFaces() + 1 * primitive.getNumNeighborCells();
   }
   /// memory is allocated on the ghost layer for each orientation (X,Y,Z,XY,XZ,XY,XYZ) and each cell
   /// most of the direction exists (num_microedges_per_edge - 1) times
   /// for the YZ orientation it is: num_microedges_per_edge
   /// for the X orientation num_microedges_per_edge - 2
   return levelinfo::num_microedges_per_edge( level ) +
          primitive.getNumNeighborFaces() * ( 3 * ( levelinfo::num_microedges_per_edge( level ) ) - 1 ) +
          primitive.getNumNeighborCells() * ( 7 * levelinfo::num_microedges_per_edge( level ) - 7 );
}

uint_t edgedof::edgeDoFMacroFaceFunctionMemorySize( const uint_t& level, const Primitive& primitive )
{
   WALBERLA_UNUSED( primitive );
   ///"inner/own" points on the face
   uint_t innerDofs =
       3 * ( ( ( levelinfo::num_microedges_per_edge( level ) + 1 ) * levelinfo::num_microedges_per_edge( level ) ) / 2 );

   ///ghost points on one adjacent tet
   uint_t GhostDoFsOneSide = 0;
   if ( primitive.getNumNeighborCells() != 0 )
   {
      /// points in the "white up" tets
      GhostDoFsOneSide += 3 * levelinfo::num_microvertices_per_face_from_width( levelinfo::num_microedges_per_edge( level ) );
      /// points from the xyz edge
      GhostDoFsOneSide += levelinfo::num_microvertices_per_face_from_width( levelinfo::num_microedges_per_edge( level ) - 1 );
      /// points on the parallel face inside the tet
      GhostDoFsOneSide += 3 * levelinfo::num_microvertices_per_face_from_width( levelinfo::num_microedges_per_edge( level ) - 1 );
   }

   return innerDofs + primitive.getNumNeighborCells() * GhostDoFsOneSide;
}

uint_t edgedof::edgeDoFMacroCellFunctionMemorySize( const uint_t& level, const Primitive& primitive )
{
   WALBERLA_UNUSED( primitive );
   return 6 * ( levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) ) ) +
          ( levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) - 1 ) );
}

unsigned long long edgedof::edgeDoFLocalFunctionMemorySize( const uint_t & level, const std::shared_ptr< PrimitiveStorage > & storage )
{
   unsigned long long mem = 0;

   for ( const auto & it : storage->getVertices() )
   {
      mem += edgeDoFMacroVertexFunctionMemorySize( level, *it.second );
   }

   for ( const auto & it : storage->getEdges() )
   {
      mem += edgeDoFMacroEdgeFunctionMemorySize( level, *it.second );
   }

   for ( const auto & it : storage->getFaces() )
   {
      mem += edgeDoFMacroFaceFunctionMemorySize( level, *it.second );
   }

   for ( const auto & it : storage->getCells() )
   {
      mem += edgeDoFMacroCellFunctionMemorySize( level, *it.second );
   }

   return mem;
}

unsigned long long edgedof::edgeDoFGlobalFunctionMemorySize( const uint_t & level, const std::shared_ptr< PrimitiveStorage > & storage )
{
   const auto memLocal = edgeDoFLocalFunctionMemorySize( level, storage );
   const auto memGlobal = walberla::mpi::allReduce( memLocal, walberla::mpi::SUM );
   return memGlobal;
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::toVector( const EdgeDoFFunction< idx_t >&       numerator,
                                             const std::shared_ptr< VectorProxy >& vec,
                                             uint_t                                level,
                                             DoFType                               flag ) const
{
   if constexpr ( !std::is_same< ValueType, real_t >::value )
   {
      WALBERLA_UNUSED( numerator );
      WALBERLA_UNUSED( vec );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( flag );
      WALBERLA_ABORT( "EdgeDoFFunction< T >::toVector() not implemented for T = " << typeid( ValueType ).name() );
   }
   else
   {
      for ( auto& it : this->getStorage()->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            edgedof::macroedge::createVectorFromFunction< real_t >(
                level, edge, this->getEdgeDataID(), numerator.getEdgeDataID(), vec );
         }
      }

      for ( auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            edgedof::macroface::createVectorFromFunction< real_t >(
                level, face, this->getFaceDataID(), numerator.getFaceDataID(), vec );
         }
      }

      for ( auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         const DoFType cellBC = this->getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
         if ( testFlag( cellBC, flag ) )
         {
            edgedof::macrocell::createVectorFromFunction< real_t >(
                level, cell, this->getCellDataID(), numerator.getCellDataID(), vec );
         }
      }
   }
}

template < typename ValueType >
void EdgeDoFFunction< ValueType >::fromVector( const EdgeDoFFunction< idx_t >&       numerator,
                                               const std::shared_ptr< VectorProxy >& vec,
                                               uint_t                                level,
                                               DoFType                               flag ) const
{
   if constexpr ( !std::is_same< ValueType, real_t >::value )
   {
      WALBERLA_UNUSED( numerator );
      WALBERLA_UNUSED( vec );
      WALBERLA_UNUSED( level );
      WALBERLA_UNUSED( flag );
      WALBERLA_ABORT( "EdgeDoFFunction< T >::fromVector() not implemented for T = " << typeid( ValueType ).name() );
   }
   else
   {
      this->startCommunication< Vertex, Edge >( level );
      this->endCommunication< Vertex, Edge >( level );

      for ( auto& it : this->getStorage()->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            edgedof::macroedge::createFunctionFromVector< real_t >(
                level, edge, this->getEdgeDataID(), numerator.getEdgeDataID(), vec );
         }
      }

      this->startCommunication< Edge, Face >( level );
      this->endCommunication< Edge, Face >( level );

      for ( auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            edgedof::macroface::createFunctionFromVector< real_t >(
                level, face, this->getFaceDataID(), numerator.getFaceDataID(), vec );
         }
      }

      for ( auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         const DoFType cellBC = this->getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
         if ( testFlag( cellBC, flag ) )
         {
            edgedof::macrocell::createFunctionFromVector< real_t >(
                level, cell, this->getCellDataID(), numerator.getCellDataID(), vec );
         }
      }
   }
}

// =================
//  specialisations
// =================
template <>
void EdgeDoFFunction< real_t >::invertElementwise( uint_t level, DoFType flag, bool workOnHalos ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Invert elementwise" );

   if ( workOnHalos )
   {
      for ( const auto& it : this->getStorage()->getVertices() )
      {
         Vertex& vertex = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
         {
            real_t* data = vertex.getData( vertexDataID_ )->getPointer( level );
            uint_t  size = edgedof::edgeDoFMacroVertexFunctionMemorySize( level, vertex );
            for ( uint_t k = 0; k < size; ++k )
            {
               data[k] = real_c( 1.0 ) / data[k];
            }
         }
      }
      for ( const auto& it : this->getStorage()->getEdges() )
      {
         Edge& edge = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
         {
            real_t* data = edge.getData( edgeDataID_ )->getPointer( level );
            uint_t  size = edgedof::edgeDoFMacroEdgeFunctionMemorySize( level, edge );
            for ( uint_t k = 0; k < size; ++k )
            {
               data[k] = real_c( 1.0 ) / data[k];
            }
         }
      }

      for ( const auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
         {
            real_t* data = face.getData( faceDataID_ )->getPointer( level );
            uint_t  size = edgedof::edgeDoFMacroFaceFunctionMemorySize( level, face );
            for ( uint_t k = 0; k < size; ++k )
            {
               data[k] = real_c( 1.0 ) / data[k];
            }
         }
      }

      for ( const auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            real_t* data = cell.getData( cellDataID_ )->getPointer( level );
            uint_t  size = edgedof::edgeDoFMacroCellFunctionMemorySize( level, cell );
            for ( uint_t k = 0; k < size; ++k )
            {
               data[k] = real_c( 1.0 ) / data[k];
            }
         }
      }
   }

   // do not work on halos
   else
   {
      for ( const auto& it : this->getStorage()->getEdges() )
      {
         Edge& edge = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
         {
            edgedof::macroedge::invertElementwise( level, edge, edgeDataID_ );
         }
      }

      for ( const auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
         {
            edgedof::macroface::invertElementwise( level, face, faceDataID_ );
         }
      }

      for ( const auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         if ( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
         {
            edgedof::macrocell::invertElementwise( level, cell, cellDataID_ );
         }
      }
   }

   this->stopTiming( "Invert elementwise" );
}

// ========================
//  explicit instantiation
// ========================
template class EdgeDoFFunction< double >;
// template class EdgeDoFFunction< float >;
template class EdgeDoFFunction< int32_t >;
template class EdgeDoFFunction< int64_t >;

template void EdgeDoFFunction< double >::interpolateByPrimitiveType< hyteg::Vertex >( const double& constant,
                                                                                      uint_t        level,
                                                                                      DoFType       flag ) const;

template void EdgeDoFFunction< double >::interpolateByPrimitiveType< hyteg::Edge >( const double& constant,
                                                                                    uint_t        level,
                                                                                    DoFType       flag ) const;

template void EdgeDoFFunction< double >::interpolateByPrimitiveType< hyteg::Face >( const double& constant,
                                                                                    uint_t        level,
                                                                                    DoFType       flag ) const;

template void EdgeDoFFunction< double >::interpolateByPrimitiveType< hyteg::Cell >( const double& constant,
                                                                                    uint_t        level,
                                                                                    DoFType       flag ) const;

} // namespace hyteg
