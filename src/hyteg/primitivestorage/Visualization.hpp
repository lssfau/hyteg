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

#include "core/Format.hpp"
#include "core/debug/CheckFunctions.h"
#include "core/mpi/MPITextFile.h"

#include "hyteg/primitives/Vertex.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

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
      rankOut << "  <VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">\n";
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
         WALBERLA_ASSERT_EQUAL( primitiveID.getID(), vertex->getID().getID() );
         vertexPosition[primitiveID.getID()] = counter;
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
            vertexPosition[vertex->getID().getID()] = counter;
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
         rankOut << vertexPosition.at( primitiveID.getID() ) << " ";
      }
      else
      {
         auto primitive = storage.getPrimitive( primitiveID );
         for ( const auto& neighborVertexID : primitive->neighborVertices() )
         {
            WALBERLA_ASSERT_GREATER( vertexPosition.count( neighborVertexID.getID() ), 0 );
            rankOut << vertexPosition.at( neighborVertexID.getID() ) << " ";
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

   rankOut << "          <DataArray type=\"UInt32\" Name=\"primitiveID\">\n";
   for ( uint_t primitive = 0; primitive < numLocalPrimitives; primitive++ )
   {
      rankOut << "            " << uint_c( primitiveIDs[primitive].getID() ) << "\n";
   }
   rankOut << "          </DataArray>\n";

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

} // namespace hyteg
