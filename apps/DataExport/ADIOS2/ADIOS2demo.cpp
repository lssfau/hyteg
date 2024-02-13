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

// This demo app allows exporting a 2D base mesh in ADIOS2 BP4 format.
//
// It works sequentially, but also in parallel. Via an XML file we can
// dynamically set e.g. the number of aggregators.
//
// NOTE: This demo was used as part of the development of HyTeG's
//       AdiosWriter. It is _not_ a tutorial. If you want to export
//       FE functions for visualisation, please use the AdiosWriter
//       class.

#include <adios2.h>
#include <sstream>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/logging/Logging.h"

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
// #include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using namespace hyteg;

void usage()
{
   WALBERLA_LOG_INFO_ON_ROOT( "Please provide the name of a 2D mesh file!" );
   WALBERLA_ABORT( "Ich kann so nicht arbeiten!" );
}

// determine number of nodes in local mesh
uint_t getNumLocalNodes( const std::shared_ptr< PrimitiveStorage > storage )
{
   std::set< PrimitiveID > nodeIDs;
   for ( const auto& it : storage->getFaces() )
   {
      const Face& face = *it.second;
      nodeIDs.insert( face.getVertexID0() );
      nodeIDs.insert( face.getVertexID1() );
      nodeIDs.insert( face.getVertexID2() );
   }
   return nodeIDs.size();
}

void addSampleData( adios2::IO& ioObject, adios2::Engine& bpWriter, uint_t numElements, uint_t numVertices )
{
   // REMARK:
   //
   // Currently, CellData seem to not be imported by ParaView (cf. MFEM example BP files). But at least they get
   // silently ignored without producing a segfault in ParaView!
   //
   // add MPI rank of this process to elements as data value
   adios2::Variable< real_t >       varMpiRank = ioObject.DefineVariable< real_t >( "MPIrank", {}, {}, { numElements } );
   adios2::Variable< real_t >::Span mpiRank    = bpWriter.Put( varMpiRank );
   for ( uint_t k = 0; k < numElements; ++k )
   {
      mpiRank[k] = static_cast< real_t >( walberla::MPIManager::instance()->rank() );
   }

   // add index of node as data values
   adios2::Variable< real_t >       varNodeIndex = ioObject.DefineVariable< real_t >( "nodeIndex", {}, {}, { numVertices } );
   adios2::Variable< real_t >::Span nodeIndex    = bpWriter.Put( varNodeIndex );
   for ( uint_t k = 0; k < numVertices; ++k )
   {
      nodeIndex[k] = static_cast< real_t >( k );
   }
}

void sequentialMeshData( const std::shared_ptr< PrimitiveStorage > storage, adios2::IO& ioObject, adios2::Engine& bpWriter )
{
   const uint_t dim         = 2;
   const uint_t numElements = storage->getNumberOfGlobalFaces();
   const uint_t numVertices = storage->getNumberOfGlobalVertices();

   // set basic entity counts
   adios2::Variable< uint32_t > varNumberOfNodes =
       ioObject.DefineVariable< uint32_t >( "NumberOfVertices", { adios2::LocalValueDim } );
   adios2::Variable< uint32_t > varNumberOfFaces =
       ioObject.DefineVariable< uint32_t >( "NumberOfElements", { adios2::LocalValueDim } );
   bpWriter.Put( varNumberOfNodes, static_cast< uint32_t >( numVertices ) );
   bpWriter.Put( varNumberOfFaces, static_cast< uint32_t >( numElements ) );

   // set element types (all the same)
   adios2::Variable< uint32_t > varElementTypes = ioObject.DefineVariable< uint32_t >( "types" );
   bpWriter.Put( varElementTypes, 5u );

   // set vertex coordinates
   adios2::Variable< real_t >       varVertices = ioObject.DefineVariable< real_t >( "vertices", {}, {}, { numVertices, dim } );
   adios2::Variable< real_t >::Span vertices    = bpWriter.Put< real_t >( varVertices );

   std::map< PrimitiveID, uint32_t > vertexIndex;
   uint_t                            idx = 0;

   for ( const auto& it : storage->getVertices() )
   {
      Point3D coords        = ( *it.second ).getCoordinates();
      vertices[2 * idx + 0] = coords( 0 );
      vertices[2 * idx + 1] = coords( 1 );
      vertexIndex[it.first] = static_cast< uint32_t >( idx );
      idx++;
   }

   // set connectivity
   adios2::Variable< uint64_t > varConnectivity =
       ioObject.DefineVariable< uint64_t >( "connectivity", {}, {}, { numElements, dim + 2 } );
   adios2::Variable< uint64_t >::Span connectivity = bpWriter.Put< uint64_t >( varConnectivity );

   idx = 0;
   for ( const auto& it : storage->getFaces() )
   {
      const Face& face      = *it.second;
      connectivity[idx + 0] = 3; // number of vertices of this element (redundant, as we also specify the type???)
      connectivity[idx + 1] = vertexIndex[face.getVertexID0()];
      connectivity[idx + 2] = vertexIndex[face.getVertexID1()];
      connectivity[idx + 3] = vertexIndex[face.getVertexID2()];
      idx += 4;
   }

   addSampleData( ioObject, bpWriter, numElements, numVertices );
}

void parallelMeshData( const std::shared_ptr< PrimitiveStorage > storage, adios2::IO& ioObject, adios2::Engine& bpWriter )
{
   const uint_t dim         = 2;
   const uint_t numElements = storage->getNumberOfLocalFaces();
   const uint_t numVertices = getNumLocalNodes( storage );

   WALBERLA_LOG_INFO( "Local mesh has " << numElements << " elements" );
   WALBERLA_LOG_INFO( "Local mesh has " << numVertices << " nodes" );

   // set basic entity counts
   adios2::Variable< uint32_t > varNumberOfNodes =
       ioObject.DefineVariable< uint32_t >( "NumberOfVertices", { adios2::LocalValueDim } );
   adios2::Variable< uint32_t > varNumberOfFaces =
       ioObject.DefineVariable< uint32_t >( "NumberOfElements", { adios2::LocalValueDim } );
   bpWriter.Put( varNumberOfNodes, static_cast< uint32_t >( numVertices ) );
   bpWriter.Put( varNumberOfFaces, static_cast< uint32_t >( numElements ) );

   // set element types (all the same)
   adios2::Variable< uint32_t > varElementTypes = ioObject.DefineVariable< uint32_t >( "types" );
   bpWriter.Put( varElementTypes, 5u );

   // set vertex coordinates and connectivities
   adios2::Variable< real_t >   varVertices = ioObject.DefineVariable< real_t >( "vertices", {}, {}, { numVertices, dim } );
   adios2::Variable< uint64_t > varConnectivity =
       ioObject.DefineVariable< uint64_t >( "connectivity", {}, {}, { numElements, dim + 2 } );

   adios2::Variable< real_t >::Span   vertices     = bpWriter.Put< real_t >( varVertices );
   adios2::Variable< uint64_t >::Span connectivity = bpWriter.Put< uint64_t >( varConnectivity );

   uint_t                          nodeIdx = 0;
   uint_t                          elemIdx = 0;
   std::map< PrimitiveID, uint_t > nodeMap;
   for ( const auto& it : storage->getFaces() )
   {
      const Face& face = *it.second;

      // number of vertices of this element (always three for triangle)
      connectivity[4 * elemIdx + 0] = 3u;

      // get PrimitiveIDs for this face's vertices
      std::vector< PrimitiveID > elemNodeIDs;
      face.getNeighborVertices( elemNodeIDs );
      WALBERLA_ASSERT( elemNodeIDs.size() == 3 );

      for ( uint_t k = 0; k < 3; ++k )
      {
         // check, if coordinates for this vertex were already stored, if not, do it
         if ( nodeMap.count( elemNodeIDs[k] ) == 0 )
         {
            const auto coords         = storage->getVertex( elemNodeIDs[k] )->getCoordinates();
            vertices[2 * nodeIdx + 0] = coords( 0 );
            vertices[2 * nodeIdx + 1] = coords( 1 );

            nodeMap[elemNodeIDs[k]] = nodeIdx;
            nodeIdx++;
         }

         // fetch node index and insert it into connectivity list
         connectivity[4 * elemIdx + 1 + k] = static_cast< uint64_t >( nodeMap[elemNodeIDs[k]] );
      }

      // proceed to next element
      elemIdx++;
   }
   WALBERLA_ASSERT( nodeIdx == numVertices );

   // add some data to visualise
   addSampleData( ioObject, bpWriter, numElements, numVertices );
}

int main( int argc, char* argv[] )
{
   // ---------------
   //  General Setup
   // ---------------
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   std::string meshFileName;

   if ( argc != 2 )
   {
      usage();
   }
   else
   {
      meshFileName = argv[1];
      WALBERLA_LOG_INFO_ON_ROOT( "Running with meshFile '" << meshFileName << "'" );
   }

   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // --------------
   //  ADIOS2 Setup
   // --------------
   std::string filename{ "output/mesh.bp" };

#ifdef WALBERLA_BUILD_WITH_MPI
   adios2::ADIOS adios( "ADIOS2config.xml", walberla::MPIManager::instance()->comm() );
#else
   adios2::ADIOS adios( "ADIOS2config.xml" );
#endif

   adios2::IO ioObject = adios.DeclareIO( "BPWriter" );
   // ioObject.SetEngine( "BP5" ); --> unsupported by ParaView 5.11.1 ?
   // ioObject.SetEngine( "BP4" ); --> set via XML file
   adios2::Engine bpWriter = ioObject.Open( filename, adios2::Mode::Write );

   // -----------------------
   //  Prepare VTU Mesh Info
   // -----------------------
   std::stringstream oStream;

   std::string myIndent1{ "  " };
   std::string myIndent2 = myIndent1 + myIndent1;
   std::string myIndent3 = myIndent1 + myIndent2;
   std::string myIndent4 = myIndent1 + myIndent3;

   std::string byteOrder = "LittleEndian"; // <-- see issue #217

   // REMARK:
   //
   // Although we specify the names of the DataArrays in the vtk.xml part, these are
   // not free to choose, but must be
   // - vertices
   // - connectivity
   // - types
   oStream << R"(<?xml version="1.0"?>)" << '\n'
           << R"(<VTKFile type="UnstructuredGrid" version="0.1" byte_order=")" << byteOrder << R"(">)" << '\n'
           << myIndent1 << R"(<UnstructuredGrid>)" << '\n'
           << myIndent2 << R"(<Piece NumberOfPoints="NumberOfVertices" NumberOfCells="NumberOfElements">)" << '\n'
           << myIndent3 << R"(<Points>)" << '\n'
           << myIndent4 << R"(<DataArray Name="vertices"/>)" << '\n'
           << myIndent3 << R"(</Points>)" << '\n'
           << myIndent3 << R"(<Cells>)" << '\n'
           << myIndent4 << R"(<DataArray Name="connectivity"/>)" << '\n'
           << myIndent4 << R"(<DataArray Name="types"/>)" << '\n'
           << myIndent3 << R"(</Cells>)" << '\n'
           << myIndent3 << R"(<PointData>)" << '\n'
           << myIndent4 << R"(<DataArray Name="nodeIndex"/>)" << '\n'
           << myIndent3 << R"(</PointData>)" << '\n'
           << myIndent3 << R"(<CellData Scalars="MPIrank">)" << '\n'
           << myIndent4 << R"(<DataArray Name="MPIrank"/>)" << '\n'
           << myIndent3 << R"(</CellData>)" << '\n'
           << myIndent2 << R"(</Piece>)" << '\n'
           << myIndent1 << R"(</UnstructuredGrid>)" << '\n'
           << R"(</VTKFile>)";

   WALBERLA_LOG_INFO_ON_ROOT( "" << oStream.rdbuf() );

   // define an attribute to store the VTU file info
   std::string attrName{ "VTU Mesh" };
   // adios2::Attribute< std::string > vtuFile = ioObject.DefineAttribute( "vtu.xml", oStream.str(), "", "" );
   ioObject.DefineAttribute( "vtk.xml", oStream.str(), "", "" );

   // -------------------
   //  Collect Mesh Data
   // -------------------
   if ( walberla::mpi::MPIManager::instance()->numProcesses() == 1 )
   {
      sequentialMeshData( storage, ioObject, bpWriter );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Detected Parallel Setting!" );
      if( storage->getNumberOfLocalFaces() > 0 ) {
        parallelMeshData( storage, ioObject, bpWriter );
      }
   }

   // ------------
   //  Perform IO
   // ------------
   bpWriter.Close();
}
