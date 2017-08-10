
#pragma once

#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "core/mpi/MPITextFile.h"

#include <memory>

namespace hhg {

void writeDomainPartitioningVTK( const std::shared_ptr< PrimitiveStorage > & storage,
                                 const std::string & dir, const std::string & filename )
{
  uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );
  uint_t numberOfProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

  auto getFilenameOfRank = []( const std::string & filename, const uint_t & rank ) -> std::string
  {
    return fmt::format("{}-rank-{:0>4}.vtu", filename, rank);
  };

  WALBERLA_ROOT_SECTION()
  {
    std::string pvtu_filename(fmt::format("{}/{}.pvtu", dir, filename));
    std::ofstream pvtu_file;
    pvtu_file.open(pvtu_filename.c_str());

    WALBERLA_CHECK( !!pvtu_file, "Error opening file: " << pvtu_filename );

    pvtu_file << "<?xml version=\"1.0\"?>\n";
    pvtu_file << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">\n";
    pvtu_file << "  <PUnstructuredGrid GhostLevel=\"0\">\n";

    // parallel point coordinates
    pvtu_file << "    <PPoints>\n";
    pvtu_file << "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\"/>\n";
    pvtu_file << "    </PPoints>\n";

    // parallel cell connectivity and types
    pvtu_file << "    <PCells>\n";
    pvtu_file << "      <PDataArray type=\"Int32\" Name=\"connectivity\"/>\n";
    pvtu_file << "      <PDataArray type=\"Int32\" Name=\"offsets\"/>\n";
    pvtu_file << "      <PDataArray type=\"UInt8\" Name=\"types\"/>\n";
    pvtu_file << "    </PCells>\n";

    // parallel cell data
    pvtu_file << "    <PCellData>\n";
    pvtu_file << "      <PDataArray type=\"UInt32\" Name=\"rank\"/>\n";
    pvtu_file << "    </PCellData>\n";

    // assemble files
    for (uint_t rankIdx = 0; rankIdx < numberOfProcesses; ++rankIdx)
    {
      pvtu_file << "   <Piece Source=\"" << getFilenameOfRank( filename, rankIdx ) << "\"/>\n";
    }

    pvtu_file << "  </PUnstructuredGrid>\n";
    pvtu_file << "</VTKFile>\n";

    pvtu_file.close();
  }

  std::string vtu_filename(fmt::format("{}/{}", dir, getFilenameOfRank( filename, rank )));
  std::ofstream vtu_file;
  vtu_file.open(vtu_filename.c_str());

  WALBERLA_CHECK( !!vtu_file, "Error opening file: " << vtu_filename );

  // header
  vtu_file << "<?xml version=\"1.0\"?>\n";
  vtu_file << "  <VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">\n";
  vtu_file << "    <UnstructuredGrid>\n";

  vtu_file << "      <Piece"
           << " NumberOfPoints=\"" << storage->getNumberOfLocalEdges() * 2 << "\""
           << " NumberOfCells=\"" << storage->getNumberOfLocalEdges() << "\""
           << ">\n";

  /////////////////////////////
  // local point coordinates //
  /////////////////////////////

  // map that maps the ID of the macro vertex to the position in the list in the .vtu file
  std::map< PrimitiveID::IDType, uint_t > vertexPosition;
  vtu_file << "        <Points>\n";
  vtu_file << "          <DataArray type=\"Float32\" NumberOfComponents=\"3\">\n";
  // write coordinates
  uint_t counter = 0;
  for ( const auto & it : storage->getEdges() )
  {
    auto edge = it.second;
    for ( const auto & neighborVertexID : edge->neighborVertices() )
    {
      WALBERLA_ASSERT(    storage->vertexExistsLocally( neighborVertexID )
                       || storage->vertexExistsInNeighborhood( neighborVertexID ));
      auto vertex = storage->getVertex( neighborVertexID );
      Point3D coordinates = vertex->getCoordinates();
      vtu_file << "            " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << "\n";
      vertexPosition[ vertex->getID().getID() ] = counter;
      counter++;
    }
  }
  vtu_file << "          </DataArray>\n";
  vtu_file << "        </Points>\n";

  ///////////////////////////////////////
  // local cell connectivity and types //
  ///////////////////////////////////////

  uint_t cellType = 3; // Edge TODO: general function for all types
  uint_t offset   = 2;

  vtu_file << "        <Cells>\n";

  // write edge connectivity
  vtu_file << "          <DataArray type=\"Int32\" Name=\"connectivity\">\n";
  for ( const auto & it : storage->getEdges() )
  {
    auto edge = it.second;
    vtu_file << "            ";
    for ( const auto & neighborVertexID : edge->neighborVertices() )
    {
      WALBERLA_ASSERT_GREATER( vertexPosition.count( neighborVertexID.getID() ), 0 );
      vtu_file << vertexPosition.at( neighborVertexID.getID() ) << " ";
    }
    vtu_file << "\n";
  }
  vtu_file << "          </DataArray>\n";

  // write offsets
  vtu_file << "          <DataArray type=\"Int32\" Name=\"offsets\">\n";
  for ( uint_t primitive = 1; primitive <= storage->getNumberOfLocalEdges(); primitive++ )
  {
    vtu_file << "            " << offset * primitive << "\n";
  }
  vtu_file << "          </DataArray>\n";

  // write cell type
  vtu_file << "          <DataArray type=\"UInt8\" Name=\"types\">\n";
  for ( uint_t primitive = 1; primitive <= storage->getNumberOfLocalEdges(); primitive++ )
  {
    vtu_file << "            " << cellType << "\n";
  }
  vtu_file << "          </DataArray>\n";

  vtu_file << "        </Cells>\n";

  /////////////////////
  // local cell data //
  /////////////////////

  vtu_file << "        <CellData>\n";
  vtu_file << "          <DataArray type=\"UInt32\" Name=\"rank\">\n";
  for ( uint_t primitive = 1; primitive <= storage->getNumberOfLocalEdges(); primitive++ )
  {
    vtu_file << "            " << rank << "\n";
  }
  vtu_file << "          </DataArray>\n";
  vtu_file << "        </CellData>\n";

  vtu_file << "      </Piece>\n";

  // closing header
  vtu_file << "    </UnstructuredGrid>\n";
  vtu_file << "  </VTKFile>\n";

  vtu_file.close();
}

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
