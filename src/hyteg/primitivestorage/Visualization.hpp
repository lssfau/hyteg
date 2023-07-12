/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include <memory>

#include "core/DataTypes.h"
#include "core/Format.hpp"
#include "core/debug/CheckFunctions.h"
#include "core/mpi/MPITextFile.h"

#include "hyteg/dataexport/VTKOutput/VTKHelpers.hpp"
#include "hyteg/primitives/Edge.hpp"
#include "hyteg/primitives/Vertex.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

enum VTK_CELL_TYPE
{
   VTK_VERTEX   = 1,
   VTK_LINE     = 3,
   VTK_TRIANGLE = 5,
   VTK_TETRA    = 10
};

static void writeDomainPartitioningVTK( const PrimitiveStorage&                                         storage,
                                        const std::string&                                              dir,
                                        const std::string&                                              filename,
                                        const VTK_CELL_TYPE&                                            vtkCellType,
                                        const std::map< std::string, std::map< PrimitiveID, real_t > >& realData )
{
   const uint_t rank              = uint_c( walberla::mpi::MPIManager::instance()->rank() );
   const uint_t numberOfProcesses = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );

   std::vector< PrimitiveID > primitiveIDs;
   uint_t                     numLocalPrimitives   = 0;
   uint_t                     verticesPerPrimitive = 0;

   switch ( vtkCellType )
   {
   case VTK_VERTEX:
      storage.getVertexIDs( primitiveIDs );
      numLocalPrimitives   = storage.getNumberOfLocalVertices();
      verticesPerPrimitive = 1;
      break;
   case VTK_LINE:
      storage.getEdgeIDs( primitiveIDs );
      numLocalPrimitives   = storage.getNumberOfLocalEdges();
      verticesPerPrimitive = 2;
      break;
   case VTK_TRIANGLE:
      storage.getFaceIDs( primitiveIDs );
      numLocalPrimitives   = storage.getNumberOfLocalFaces();
      verticesPerPrimitive = 3;
      break;
   case VTK_TETRA:
      storage.getCellIDs( primitiveIDs );
      numLocalPrimitives   = storage.getNumberOfLocalCells();
      verticesPerPrimitive = 4;
      break;
   default:
      WALBERLA_ASSERT( false, "VTK cell type not supported!" );
      break;
   }

   std::ostringstream rankOut;

   WALBERLA_ROOT_SECTION()
   {
      // header
      rankOut << "<?xml version=\"1.0\"?>\n";
      rankOut << "  <VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"" << vtk::getByteOrder() << "\">\n";
      rankOut << "    <UnstructuredGrid>\n";
   }

   rankOut << "      <Piece"
           << " NumberOfPoints=\"" << numLocalPrimitives * verticesPerPrimitive << "\""
           << " NumberOfCells=\"" << numLocalPrimitives << "\""
           << ">\n";

   /////////////////////////////
   // local point coordinates //
   /////////////////////////////

   // map that maps the ID of the macro vertex to the position in the list in the .vtu file
   std::map< PrimitiveID, uint_t > vertexPosition;
   rankOut << "        <Points>\n";
   rankOut << "          <DataArray type=\"Float32\" NumberOfComponents=\"3\">\n";
   // write coordinates
   uint_t counter = 0;
   for ( const auto& primitiveID : primitiveIDs )
   {
      if ( vtkCellType == VTK_VERTEX )
      {
         auto    vertex      = storage.getVertex( primitiveID );
         Point3D coordinates = vertex->getCoordinates();
         rankOut << "            " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << "\n";
         WALBERLA_ASSERT_EQUAL( primitiveID, vertex->getID() );
         vertexPosition[primitiveID] = counter;
         counter++;
      }
      else
      {
         auto primitive = storage.getPrimitive( primitiveID );
         WALBERLA_ASSERT_EQUAL( primitive->getNumNeighborVertices(), verticesPerPrimitive );
         for ( const auto& neighborVertexID : primitive->neighborVertices() )
         {
            WALBERLA_ASSERT( storage.vertexExistsLocally( neighborVertexID ) ||
                             storage.vertexExistsInNeighborhood( neighborVertexID ) );
            auto    vertex      = storage.getVertex( neighborVertexID );
            Point3D coordinates = vertex->getCoordinates();
            rankOut << "            " << coordinates[0] << " " << coordinates[1] << " " << coordinates[2] << "\n";
            vertexPosition[vertex->getID()] = counter;
            counter++;
         }
      }
   }
   rankOut << "          </DataArray>\n";
   rankOut << "        </Points>\n";

   ///////////////////////////////////////
   // local cell connectivity and types //
   ///////////////////////////////////////

   const uint_t offset = verticesPerPrimitive;

   rankOut << "        <Cells>\n";

   // write connectivity
   rankOut << "          <DataArray type=\"Int32\" Name=\"connectivity\">\n";
   for ( const auto& primitiveID : primitiveIDs )
   {
      rankOut << "            ";
      if ( vtkCellType == VTK_VERTEX )
      {
         rankOut << vertexPosition.at( primitiveID ) << " ";
      }
      else
      {
         auto primitive = storage.getPrimitive( primitiveID );
         for ( const auto& neighborVertexID : primitive->neighborVertices() )
         {
            WALBERLA_ASSERT_GREATER( vertexPosition.count( neighborVertexID ), 0 );
            rankOut << vertexPosition.at( neighborVertexID ) << " ";
         }
      }
      rankOut << "\n";
   }
   rankOut << "          </DataArray>\n";

   // write offsets
   rankOut << "          <DataArray type=\"Int32\" Name=\"offsets\">\n";
   for ( uint_t primitive = 1; primitive <= numLocalPrimitives; primitive++ )
   {
      rankOut << "            " << offset * primitive << "\n";
   }
   rankOut << "          </DataArray>\n";

   // write cell type
   rankOut << "          <DataArray type=\"UInt8\" Name=\"types\">\n";
   for ( uint_t primitive = 1; primitive <= numLocalPrimitives; primitive++ )
   {
      rankOut << "            " << (uint_t) vtkCellType << "\n";
   }
   rankOut << "          </DataArray>\n";

   rankOut << "        </Cells>\n";

   /////////////////////
   // local cell data //
   /////////////////////

   rankOut << "        <CellData>\n";
   rankOut << "          <DataArray type=\"UInt32\" Name=\"rank\">\n";
   for ( uint_t primitive = 0; primitive < numLocalPrimitives; primitive++ )
   {
      rankOut << "            " << rank << "\n";
   }
   rankOut << "          </DataArray>\n";

   rankOut << "          <DataArray type=\"UInt32\" Name=\"meshBoundaryFlag\">\n";
   for ( uint_t primitive = 0; primitive < numLocalPrimitives; primitive++ )
   {
      rankOut << "            " << storage.getPrimitive( primitiveIDs[primitive] )->getMeshBoundaryFlag() << "\n";
   }
   rankOut << "          </DataArray>\n";

   rankOut << "          <DataArray type=\"UInt32\" Name=\"level\">\n";
   for ( uint_t primitive = 0; primitive < numLocalPrimitives; primitive++ )
   {
      rankOut << "            " << storage.getRefinementLevel( primitiveIDs[primitive] ) << "\n";
   }
   rankOut << "          </DataArray>\n";

#if 0
   rankOut << "          <DataArray type=\"UInt32\" Name=\"primitiveID\">\n";
   for ( uint_t primitive = 0; primitive < numLocalPrimitives; primitive++ )
   {
      rankOut << "            " << uint_c( primitiveIDs[primitive].getID() ) << "\n";
   }
   rankOut << "          </DataArray>\n";
#endif

   for ( const auto& it : realData )
   {
      const auto  realDataName    = it.first;
      const auto& realDataContent = it.second;

      rankOut << "          <DataArray type=\"Float32\" Name=\"" << realDataName << "\">\n";
      for ( uint_t primitive = 0; primitive < numLocalPrimitives; primitive++ )
      {
         rankOut << "            " << realDataContent.at( primitiveIDs[primitive] ) << "\n";
      }
      rankOut << "          </DataArray>\n";
   }

   rankOut << "        </CellData>\n";

   rankOut << "      </Piece>\n";

   // write in parallel
   std::string vtu_filename( walberla::format( "%s/%s.vtu", dir.c_str(), filename.c_str() ) );
   walberla::mpi::writeMPITextFile( vtu_filename, rankOut.str() );

   WALBERLA_ROOT_SECTION()
   {
      std::ofstream pvtu_file;
      pvtu_file.open( vtu_filename.c_str(), std::ofstream::out | std::ofstream::app );
      WALBERLA_CHECK( !!pvtu_file, "[VTKWriter (domain partitioning)] Error opening file: " << vtu_filename );
      pvtu_file << "    </UnstructuredGrid>\n";
      pvtu_file << "  </VTKFile>\n";
      pvtu_file.close();
   }
}

inline void writeDomainPartitioningVTK( const PrimitiveStorage&                                         storage,
                                        const std::string&                                              dir,
                                        const std::string&                                              filename,
                                        const std::map< std::string, std::map< PrimitiveID, real_t > >& realData )
{
   const std::string filenameVertices = filename + "_vertices";
   const std::string filenameEdges    = filename + "_edges";
   const std::string filenameFaces    = filename + "_faces";
   const std::string filenameCells    = filename + "_cells";

   writeDomainPartitioningVTK( storage, dir, filenameVertices, VTK_VERTEX, realData );
   writeDomainPartitioningVTK( storage, dir, filenameEdges, VTK_LINE, realData );
   writeDomainPartitioningVTK( storage, dir, filenameFaces, VTK_TRIANGLE, realData );
   if ( storage.hasGlobalCells() )
   {
      writeDomainPartitioningVTK( storage, dir, filenameCells, VTK_TETRA, realData );
   }
}

inline void writeDomainPartitioningVTK( const PrimitiveStorage& storage, const std::string& dir, const std::string& filename )
{
   std::map< std::string, std::map< PrimitiveID, real_t > > realData;
   writeDomainPartitioningVTK( storage, dir, filename, realData );
}

inline void writeDomainPartitioningVTK( const std::shared_ptr< PrimitiveStorage >& storage,
                                        const std::string&                         dir,
                                        const std::string&                         filename )
{
   writeDomainPartitioningVTK( *storage, dir, filename );
}

inline void writeDomainPartitioningVTK( const std::shared_ptr< PrimitiveStorage >&                      storage,
                                        const std::string&                                              dir,
                                        const std::string&                                              filename,
                                        const std::map< std::string, std::map< PrimitiveID, real_t > >& realData )
{
   writeDomainPartitioningVTK( *storage, dir, filename, realData );
}

inline void writePrimitiveStorageDistributionCSV( const std::shared_ptr< PrimitiveStorage >& storage,
                                                  const std::string&                         filename )
{
   uint_t rank = uint_c( walberla::mpi::MPIManager::instance()->rank() );

   const std::string delimiter = " ";
   std::stringstream output;

   WALBERLA_ROOT_SECTION()
   {
      output << "processRank" << delimiter << "numPrimitives" << delimiter << "numVertices" << delimiter << "numEdges"
             << delimiter << "numFaces" << delimiter << "fractionOfLocalNeighbors" << delimiter << "\n";
   }

   uint_t numTotalLocalNeighbors = 0;
   uint_t numLocalLocalNeighbors = 0;

   PrimitiveStorage::PrimitiveMap primitiveMap;
   storage->getPrimitives( primitiveMap );

   for ( const auto& primitive : primitiveMap )
   {
      numTotalLocalNeighbors += primitive.second->getNumNeighborVertices() + primitive.second->getNumNeighborEdges() +
                                primitive.second->getNumNeighborFaces();
      std::vector< PrimitiveID > neighbors;
      primitive.second->getNeighborPrimitives( neighbors );
      for ( const auto& neighbor : neighbors )
      {
         if ( storage->primitiveExistsLocally( neighbor ) )
         {
            numLocalLocalNeighbors++;
         }
      }
   }

   real_t fractionOfLocalNeighbors = (real_t) numLocalLocalNeighbors / (real_t) numTotalLocalNeighbors;

   output << rank << delimiter << storage->getNumberOfLocalPrimitives() << delimiter << storage->getNumberOfLocalVertices()
          << delimiter << storage->getNumberOfLocalEdges() << delimiter << storage->getNumberOfLocalFaces() << delimiter
          << fractionOfLocalNeighbors << delimiter << "\n";

   walberla::mpi::writeMPITextFile( filename, output.str(), walberla::mpi::MPIManager::instance()->comm() );
}

/// \brief Writes the coarse mesh as it appears in physical space.
///
/// This function takes all macro edges and places `resolution` many vertices on them (evenly spaced in computational space).
/// It then transforms all vertices by the blending function and writes the edges to the specified VTK file.
/// Meshes on finer levels can be written by refining the mesh first with MeshInfo::refinedCoarseMesh.
///
/// \param storage    PrimitiveStorage defining the mesh.
/// \param dir        Directory where the files are stored.
/// \param filename   Basename of the VTK files.
/// \param resolution Number of vertices to put on each edge. Must be >= 2.
inline void writeBlendedCoarseMeshVTK( const PrimitiveStorage& storage,
                                       const std::string&      dir,
                                       const std::string&      filename,
                                       const uint_t            resolution )
{
   WALBERLA_ASSERT_GREATER_EQUAL( resolution, 2 );

   const uint_t numLocalEdges = storage.getNumberOfLocalEdges();

   std::ostringstream rankOut;

   WALBERLA_ROOT_SECTION()
   {
      // header
      rankOut << "<?xml version=\"1.0\"?>\n";
      rankOut << "  <VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"" << vtk::getByteOrder() << "\">\n";
      rankOut << "    <PolyData>\n";
   }

   rankOut << "      <Piece"
           << " NumberOfPoints=\"" << numLocalEdges * resolution << "\""
           << " NumberOfVerts=\"0\""
           << " NumberOfLines=\"" << numLocalEdges << "\""
           << " NumberOfStrips=\"0\""
           << " NumberOfPolys=\"0\""
           << ">\n";

   /////////////////////////////
   // local point coordinates //
   /////////////////////////////

   rankOut << "        <Points>\n";
   rankOut << "          <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n";

   // write coordinates
   for ( const auto& edgeIt : storage.getEdges() )
   {
      const std::shared_ptr< Edge > edge   = edgeIt.second;
      const auto                    coords = edge->getCoordinates();

      for ( uint_t i = 0; i < resolution; ++i )
      {
         const double  t     = walberla::numeric_cast< double >( i ) / walberla::numeric_cast< double >( resolution - 1 );
         const Point3D xComp = ( 1 - t ) * coords[0] + t * coords[1];
         Point3D       xPhys;
         edge->getGeometryMap()->evalF( xComp, xPhys );

         rankOut << "            " << xPhys[0] << " " << xPhys[1] << " " << xPhys[2] << "\n";
      }
   }
   rankOut << "          </DataArray>\n";
   rankOut << "        </Points>\n";

   /////////////////
   // local edges //
   /////////////////

   rankOut << "        <Lines>\n";

   // write connectivity
   rankOut << "          <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
   for ( uint_t edgeIdx = 0; edgeIdx < numLocalEdges; edgeIdx++ )
   {
      rankOut << "           ";
      for ( uint_t i = 0; i < resolution; i++ )
      {
         rankOut << " " << edgeIdx * resolution + i;
      }
      rankOut << "\n";
   }
   rankOut << "          </DataArray>\n";

   // write offsets
   rankOut << "          <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
   for ( uint_t i = 0; i < numLocalEdges; i++ )
   {
      rankOut << "            " << ( i + 1 ) * resolution << "\n";
   }
   rankOut << "          </DataArray>\n";

   rankOut << "        </Lines>\n";

   rankOut << "      </Piece>\n";

   // write in parallel
   std::string vtp_filename( walberla::format( "%s/%s.vtp", dir.c_str(), filename.c_str() ) );
   walberla::mpi::writeMPITextFile( vtp_filename, rankOut.str() );

   WALBERLA_ROOT_SECTION()
   {
      std::ofstream pvtp_file;
      pvtp_file.open( vtp_filename.c_str(), std::ofstream::out | std::ofstream::app );
      WALBERLA_CHECK( !!pvtp_file, "[VTKWriter (blended coarse mesh)] Error opening file: " << vtp_filename );
      pvtp_file << "    </PolyData>\n";
      pvtp_file << "  </VTKFile>\n";
      pvtp_file.close();
   }
}

/// See void writeBlendedCoarseMeshVTK( const PrimitiveStorage& storage, const std::string& dir, const std::string& filename, const uint_t resolution ).
inline void writeBlendedCoarseMeshVTK( const std::shared_ptr< PrimitiveStorage >& storage,
                                       const std::string&                         dir,
                                       const std::string&                         filename,
                                       const uint_t                               resolution = 9 )
{
   writeBlendedCoarseMeshVTK( *storage, dir, filename, resolution );
}

} // namespace hyteg
