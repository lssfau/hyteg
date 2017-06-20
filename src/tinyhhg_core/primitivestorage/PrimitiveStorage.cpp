
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
  std::map< PrimitiveID::IDType, SetupVertex* > setupVertices;
  std::map< PrimitiveID::IDType, SetupEdge* >   setupEdges;
  std::map< PrimitiveID::IDType, SetupFace* >   setupFaces;

  setupStorage.getVertices( setupVertices );
  setupStorage.getEdges( setupEdges );
  setupStorage.getFaces( setupFaces );

  uint_t processRank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

  for ( auto it = setupVertices.begin(); it != setupVertices.end(); it++  )
  {
    if ( processRank == it->second->getTargetRank() )
    {
      vertices_[ it->first ] = new Vertex( it->first, it->second->getCoordinates() );
    }
  }

  for ( auto it = setupEdges.begin(); it != setupEdges.end(); it++ )
  {
    if ( processRank == it->second->getTargetRank() )
    {
      PrimitiveID edgeID = it->first;
      DoFType edgeType = it->second->getDoFType();
      Vertex* edgeVertex0 = vertices_[ it->second->getVertexID0().getID() ];
      Vertex* edgeVertex1 = vertices_[ it->second->getVertexID1().getID() ];
      edges_[ edgeID.getID() ] = new Edge( edgeID.getID(), edgeType, edgeVertex0, edgeVertex1 );
    }
  }

  for ( auto it = setupFaces.begin(); it != setupFaces.end(); it++ )
  {
    if ( processRank == it->second->getTargetRank() )
    {
      PrimitiveID faceID = it->first;
      Edge* faceEdges[3];
      faceEdges[0] = edges_[ it->second->getEdgeID0().getID() ];
      faceEdges[1] = edges_[ it->second->getEdgeID1().getID() ];
      faceEdges[2] = edges_[ it->second->getEdgeID2().getID() ];
      faces_[ faceID.getID() ] = new Face( faceID.getID(), faceEdges );
    }
  }

#ifndef NDEBUG
  checkConsistency();
#endif
}


void PrimitiveStorage::checkConsistency()
{
  // 1. Number of data handlers equal to local counter
  // 2. PrimitiveIDs of maps match IDs of Primtives
  for ( auto it = vertices_.begin(); it != vertices_.end(); it++ )
  {
    WALBERLA_CHECK_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
  }
  for ( auto it = edges_.begin(); it != edges_.end(); it++ )
  {
    WALBERLA_CHECK_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
  }
  for ( auto it = faces_.begin(); it != faces_.end(); it++ )
  {
    WALBERLA_CHECK_EQUAL( primitiveDataHandlers_, it->second->getNumberOfDataEntries() );
    WALBERLA_CHECK_EQUAL( it->first, it->second->getID().getID() );
  }



}


} // namespace hhg

