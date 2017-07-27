
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"
#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitives/face.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"


#include <map>
#include <vector>

namespace hhg {

using walberla::uint_t;

PrimitiveStorage::PrimitiveStorage( const uint_t & rank, const SetupPrimitiveStorage & setupStorage ) :
  rank_( rank ), primitiveDataHandlers_( 0 )
{
  for ( auto it = setupStorage.beginVertices(); it != setupStorage.endVertices(); it++  )
  {
    if ( rank_ == setupStorage.getTargetRank( it->first ) )
    {
      vertices_[ it->first ] = new Vertex( *it->second );
    }
  }

  for ( auto it = setupStorage.beginEdges(); it != setupStorage.endEdges(); it++ )
  {
    if ( rank_ == setupStorage.getTargetRank( it->first ) )
    {
      edges_[ it->first ] = new Edge( *it->second );
    }
  }

  for ( auto it = setupStorage.beginFaces(); it != setupStorage.endFaces(); it++ )
  {
    if ( rank_ == setupStorage.getTargetRank( it->first ) )
    {
      faces_[ it->first ] = new Face( *it->second );
    }
  }

  // Neighborhood

  for ( const auto & it : vertices_ )
  {
    Vertex * vertex = it.second;

    for ( const auto & neighborVertexID : vertex->neighborVertices() )
    {
      const Vertex * neighborVertex = setupStorage.getVertex( neighborVertexID );
      if ( !vertexExistsLocally( neighborVertexID ) && !vertexExistsInNeighborhood( neighborVertexID ) )
      {
        neighborVertices_[ neighborVertexID.getID() ] = new Vertex( *neighborVertex );
        neighborRanks_[ neighborVertexID.getID() ] = setupStorage.getTargetRank( neighborVertexID.getID() );
      }
    }

    for ( const auto & neighborEdgeID : vertex->neighborEdges() )
    {
      const Edge * neighborEdge = setupStorage.getEdge( neighborEdgeID );
      if ( !edgeExistsLocally( neighborEdgeID ) && !edgeExistsInNeighborhood( neighborEdgeID ) )
      {
        neighborEdges_[ neighborEdgeID.getID() ] = new Edge( *neighborEdge );
        neighborRanks_[ neighborEdgeID.getID() ] = setupStorage.getTargetRank( neighborEdgeID.getID() );
      }
    }

    for ( const auto & neighborFaceID : vertex->neighborFaces() )
    {
      const Face * neighborFace = setupStorage.getFace( neighborFaceID );
      if ( !faceExistsLocally( neighborFaceID ) && !faceExistsInNeighborhood( neighborFaceID ) )
      {
        neighborFaces_[ neighborFaceID.getID() ] = new Face( *neighborFace );
        neighborRanks_[ neighborFaceID.getID() ] = setupStorage.getTargetRank( neighborFaceID.getID() );
      }
    }
  }

  for ( const auto & it : edges_ )
  {
    Edge * edge = it.second;

    for ( const auto & neighborVertexID : edge->neighborVertices() )
    {
      const Vertex * neighborVertex = setupStorage.getVertex( neighborVertexID );
      if ( !vertexExistsLocally( neighborVertexID ) && !vertexExistsInNeighborhood( neighborVertexID ) )
      {
        neighborVertices_[ neighborVertexID.getID() ] = new Vertex( *neighborVertex );
        neighborRanks_[ neighborVertexID.getID() ] = setupStorage.getTargetRank( neighborVertexID.getID() );
      }
    }

    for ( const auto & neighborEdgeID : edge->neighborEdges() )
    {
      const Edge * neighborEdge = setupStorage.getEdge( neighborEdgeID );
      if ( !edgeExistsLocally( neighborEdgeID ) && !edgeExistsInNeighborhood( neighborEdgeID ) )
      {
        neighborEdges_[ neighborEdgeID.getID() ] = new Edge( *neighborEdge );
        neighborRanks_[ neighborEdgeID.getID() ] = setupStorage.getTargetRank( neighborEdgeID.getID() );
      }
    }

    for ( const auto & neighborFaceID : edge->neighborFaces() )
    {
      const Face * neighborFace = setupStorage.getFace( neighborFaceID );
      if ( !faceExistsLocally( neighborFaceID ) && !faceExistsInNeighborhood( neighborFaceID ) )
      {
        neighborFaces_[ neighborFaceID.getID() ] = new Face( *neighborFace );
        neighborRanks_[ neighborFaceID.getID() ] = setupStorage.getTargetRank( neighborFaceID.getID() );
      }
    }
  }

  for ( const auto & it : faces_ )
  {
    Face * face = it.second;

    for ( const auto & neighborVertexID : face->neighborVertices() )
    {
      const Vertex * neighborVertex = setupStorage.getVertex( neighborVertexID );
      if ( !vertexExistsLocally( neighborVertexID ) && !vertexExistsInNeighborhood( neighborVertexID ) )
      {
        neighborVertices_[ neighborVertexID.getID() ] = new Vertex( *neighborVertex );
        neighborRanks_[ neighborVertexID.getID() ] = setupStorage.getTargetRank( neighborVertexID.getID() );
      }
    }

    for ( const auto & neighborEdgeID : face->neighborEdges() )
    {
      const Edge * neighborEdge = setupStorage.getEdge( neighborEdgeID );
      if ( !edgeExistsLocally( neighborEdgeID ) && !edgeExistsInNeighborhood( neighborEdgeID ) )
      {
        neighborEdges_[ neighborEdgeID.getID() ] = new Edge( *neighborEdge );
        neighborRanks_[ neighborEdgeID.getID() ] = setupStorage.getTargetRank( neighborEdgeID.getID() );
      }
    }

    for ( const auto & neighborFaceID : face->neighborFaces() )
    {
      const Face * neighborFace = setupStorage.getFace( neighborFaceID );
      if ( !faceExistsLocally( neighborFaceID ) && !faceExistsInNeighborhood( neighborFaceID ) )
      {
        neighborFaces_[ neighborFaceID.getID() ] = new Face( *neighborFace );
        neighborRanks_[ neighborFaceID.getID() ] = setupStorage.getTargetRank( neighborFaceID.getID() );
      }
    }
  }


#ifndef NDEBUG
  checkConsistency();
#endif
}



void PrimitiveStorage::getPrimitives( PrimitiveMap & primitiveMap ) const
{
  primitiveMap.clear();

  primitiveMap.insert( beginVertices(), endVertices() );
  primitiveMap.insert( beginEdges(), endEdges() );
  primitiveMap.insert( beginFaces(), endFaces() );

  WALBERLA_ASSERT_EQUAL( primitiveMap.size(), vertices_.size() + edges_.size() + faces_.size() );
}


const Primitive* PrimitiveStorage::getPrimitive( const PrimitiveID & id ) const
{
  WALBERLA_ASSERT_LESS_EQUAL(   vertices_.count( id.getID() )
                              +    edges_.count( id.getID() )
                              +    faces_.count( id.getID() ), uint_c( 1 ) );
  if ( vertexExistsLocally( id ) ) return getVertex( id );
  if (   edgeExistsLocally( id ) ) return   getEdge( id );
  if (   faceExistsLocally( id ) ) return   getFace( id );
  return nullptr;
}

Primitive* PrimitiveStorage::getPrimitive( const PrimitiveID & id )
{
  WALBERLA_ASSERT_LESS_EQUAL(   vertices_.count( id.getID() )
                              +    edges_.count( id.getID() )
                              +    faces_.count( id.getID() ), uint_c( 1 ) );
  if ( vertexExistsLocally( id ) ) return getVertex( id );
  if (   edgeExistsLocally( id ) ) return   getEdge( id );
  if (   faceExistsLocally( id ) ) return   getFace( id );
  return nullptr;
}


const Primitive* PrimitiveStorage::getNeighborPrimitive( const PrimitiveID & id ) const
{
  WALBERLA_ASSERT_LESS_EQUAL(   neighborVertices_.count( id.getID() )
                              +    neighborEdges_.count( id.getID() )
                              +    neighborFaces_.count( id.getID() ), uint_c( 1 ) );
  if ( vertexExistsInNeighborhood( id ) ) return getNeighborVertex( id );
  if (   edgeExistsInNeighborhood( id ) ) return   getNeighborEdge( id );
  if (   faceExistsInNeighborhood( id ) ) return   getNeighborFace( id );
  return nullptr;
}

Primitive* PrimitiveStorage::getNeighborPrimitive( const PrimitiveID & id )
{
  WALBERLA_ASSERT_LESS_EQUAL(   neighborVertices_.count( id.getID() )
                              +    neighborEdges_.count( id.getID() )
                              +    neighborFaces_.count( id.getID() ), uint_c( 1 ) );
  if ( vertexExistsInNeighborhood( id ) ) return getNeighborVertex( id );
  if (   edgeExistsInNeighborhood( id ) ) return   getNeighborEdge( id );
  if (   faceExistsInNeighborhood( id ) ) return   getNeighborFace( id );
  return nullptr;
}


void PrimitiveStorage::getVertexIDs ( std::vector< PrimitiveID > & vertexIDs ) const
{
  vertexIDs.clear();
  for ( auto const & it : vertices_ )
  {
    vertexIDs.push_back( it.first );
  }
}

void PrimitiveStorage::getEdgeIDs ( std::vector< PrimitiveID > & edgeIDs ) const
{
  edgeIDs.clear();
  for ( auto const & it : edges_ )
  {
    edgeIDs.push_back( it.first );
  }
}

void PrimitiveStorage::getFaceIDs ( std::vector< PrimitiveID > & faceIDs ) const
{
  faceIDs.clear();
  for ( auto const & it : faces_ )
  {
    faceIDs.push_back( it.first );
  }
}

uint_t PrimitiveStorage::getPrimitiveRank ( const PrimitiveID & id ) const
{
  WALBERLA_ASSERT( primitiveExistsLocally( id ) || primitiveExistsInNeighborhood( id ) );
  if ( primitiveExistsLocally( id ) )
  {
    return uint_c( walberla::mpi::MPIManager::instance()->rank() );
  }
  else
  {
    return getNeighborPrimitiveRank( id );
  }
}

void PrimitiveStorage::checkConsistency()
{
  // 1. Number of data handlers less than local counter
  // 2. PrimitiveIDs of maps match IDs of Primitives
  // 3. Neighborhood of Primitives
  for ( auto it = vertices_.begin(); it != vertices_.end(); it++ )
  {
    WALBERLA_CHECK_GREATER_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 0 );
  }
  for ( auto it = edges_.begin(); it != edges_.end(); it++ )
  {
    WALBERLA_CHECK_GREATER_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 2 );
  }
  for ( auto it = faces_.begin(); it != faces_.end(); it++ )
  {
    WALBERLA_CHECK_GREATER_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 3 );
  }

  // 4. Number of data handlers of neighbor primitives is zero
  // 5. PrimitiveIDs of neighbor maps match IDs of neighbor Primitives
  // 6. Neighborhood of Primitives
  for ( auto it = neighborVertices_.begin(); it != neighborVertices_.end(); it++ )
  {
    WALBERLA_CHECK_GREATER_EQUAL( 0, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 0 );
  }
  for ( auto it = neighborEdges_.begin(); it != neighborEdges_.end(); it++ )
  {
    WALBERLA_CHECK_GREATER_EQUAL( 0, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 2 );
  }
  for ( auto it = neighborFaces_.begin(); it != neighborFaces_.end(); it++ )
  {
    WALBERLA_CHECK_GREATER_EQUAL( 0, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
    WALBERLA_CHECK_EQUAL( it->second->getNumLowerDimNeighbors(), 3 );
  }


}


} // namespace hhg

