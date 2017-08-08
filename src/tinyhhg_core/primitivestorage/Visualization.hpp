
#pragma once

#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "core/mpi/MPITextFile.h"

#include <memory>

namespace hhg {

void writePrimitiveStorageDistributionCSV( const std::shared_ptr< PrimitiveStorage > & storage, const std::string & filename )
{
  uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

  const std::string delimiter = " ";
  std::stringstream output;

  WALBERLA_ROOT_SECTION()
  {
    output << "processRank" << delimiter
           << "numVertices" << delimiter
           << "numEdges"    << delimiter
           << "numFaces"    << delimiter
           << "fractionOfLocalVertexNeighbors" << delimiter
           << "\n";
  }

  uint_t numTotalNeighborsOfLocalVertices = 0;
  uint_t numLocalNeighborsOfLocalVertices = 0;

  for ( const auto & vertex : storage->getVertices() )
  {
    WALBERLA_ASSERT_EQUAL( vertex.second->getNumNeighborVertices(), 0 );
    numTotalNeighborsOfLocalVertices += vertex.second->getNumNeighborVertices() + vertex.second->getNumNeighborEdges() + vertex.second->getNumNeighborFaces();
    std::vector< PrimitiveID > neighbors;
    vertex.second->getNeighborPrimitives( neighbors );
    for ( const auto & neighbor : neighbors ) {
      if ( storage->primitiveExistsLocally( neighbor ) )
      {
        numLocalNeighborsOfLocalVertices++;
      }
    }
  }

  real_t fractionOfLocalVertexNeighbors = (real_t) numLocalNeighborsOfLocalVertices / (real_t) numTotalNeighborsOfLocalVertices;

  output << rank << delimiter
         << storage->getNumberOfLocalVertices() << delimiter
         << storage->getNumberOfLocalEdges()    << delimiter
         << storage->getNumberOfLocalFaces()    << delimiter
         << fractionOfLocalVertexNeighbors      << delimiter
         << "\n";

  walberla::mpi::writeMPITextFile( filename, output.str(), walberla::mpi::MPIManager::instance()->comm() );
}

}
