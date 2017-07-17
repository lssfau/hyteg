
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"

#include "core/debug/CheckFunctions.h"
#include "core/debug/Debug.h"
#include "core/logging/Logging.h"
#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"
#include "tinyhhg_core/primitives/vertex.hpp"
#include "tinyhhg_core/primitives/edge.hpp"
#include "tinyhhg_core/primitives/face.hpp"

#include <map>
#include <vector>

namespace hhg {

using walberla::uint_t;

PrimitiveStorage::PrimitiveStorage( const uint_t & rank, const SetupPrimitiveStorage & setupStorage ) :
  rank_( rank ), primitiveDataHandlers_( 0 )
{
  for ( auto it = setupStorage.beginVertices(); it != setupStorage.endVertices(); it++  )
  {
    if ( rank_ == it->second->getTargetRank() )
    {
      vertices_[ it->first ] = new Vertex( *this, setupStorage, it->first );
    }
  }

  for ( auto it = setupStorage.beginEdges(); it != setupStorage.endEdges(); it++ )
  {
    if ( rank_ == it->second->getTargetRank() )
    {
      edges_[ it->first ] = new Edge( *this, setupStorage, it->first );
    }
  }

  for ( auto it = setupStorage.beginFaces(); it != setupStorage.endFaces(); it++ )
  {
    if ( rank_ == it->second->getTargetRank() )
    {
      faces_[ it->first ] = new Face( *this, setupStorage, it->first );
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
  return NULL;
}

Primitive* PrimitiveStorage::getPrimitive( const PrimitiveID & id )
{
  WALBERLA_ASSERT_LESS_EQUAL(   vertices_.count( id.getID() )
                              +    edges_.count( id.getID() )
                              +    faces_.count( id.getID() ), uint_c( 1 ) );
  if ( vertexExistsLocally( id ) ) return getVertex( id );
  if (   edgeExistsLocally( id ) ) return   getEdge( id );
  if (   faceExistsLocally( id ) ) return   getFace( id );
  return NULL;
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



}


} // namespace hhg

