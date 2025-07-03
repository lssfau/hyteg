/*
 * Copyright (c) 2017-2024 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/dataexport/VTKOutput/VTKP2PlusBubbleWriter.hpp"

#include "core/DataTypes.h"

#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"

#include "vtk/UtilityFunctions.h"

#ifdef HYTEG_BUILD_WITH_ADIOS2
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#endif

namespace hyteg {

using walberla::vtk::typeToString;

void VTKP2PlusBubbleWriter::write( const VTKOutput& mgr, std::ostream& output, const uint_t& level )
{
   if ( mgr.getNumRegisteredFunctions( vtk::DoFType::P2_PLUS_BUBBLE ) == 0 )
   {
      return;
   }

   auto storage = mgr.storage_;

   const uint_t faceCountLevel = level;
   const uint_t cellCountLevel = level;

   const uint_t numMicroElementsPerMacro =
       mgr.write2D_ ? levelinfo::num_microfaces_per_face( faceCountLevel ) : levelinfo::num_microcells_per_cell( cellCountLevel );

   const uint_t numberOfPoints =
       mgr.write2D_ ?
           storage->getNumberOfLocalFaces() * ( levelinfo::num_microvertices_per_face( level + 1 ) + numMicroElementsPerMacro ) :
           storage->getNumberOfLocalCells() * ( levelinfo::num_microvertices_per_cell( level + 1 ) + numMicroElementsPerMacro );

   const uint_t numberOfCells = mgr.write2D_ ?
                                    storage->getNumberOfLocalFaces() * levelinfo::num_microfaces_per_face( faceCountLevel ) * 6 :
                                    storage->getNumberOfLocalCells() * levelinfo::num_microcells_per_cell( cellCountLevel ) *
                                        0; // <--- need to figure out factor

   vtk::writePieceHeader( output, numberOfPoints, numberOfCells );

   output << "<Points>\n";
   vtk::openDataElement( output, typeToString< real_t >(), "", 3, mgr.vtkDataFormat_ );

   {
      VTKStreamWriter< real_t > streamWriter( mgr.vtkDataFormat_ );

      // Locations of VertexDoFs
      VTKMeshWriter::writePointsForMicroVertices( mgr.write2D_, streamWriter, storage, level );

      // Locations of EdgeDoFs
      if ( mgr.write2D_ )
      {
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_X );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_XY );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_Y );
      }
      else
      {
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_X );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_Y );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_Z );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_XY );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_XZ );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_YZ );
         VTKMeshWriter::writePointsForMicroEdges( mgr.write2D_, streamWriter, storage, level, vtk::DoFType::EDGE_XYZ );
      }

      // Locations of Face/VolumeDoFs
      if ( mgr.write2D_ )
      {
         VTKMeshWriter::writePointsForMicroFaceCenters( streamWriter, storage, level );
      }
      else
      {
         WALBERLA_ABORT( "Missing 3D implementation inside VTKP2PlusBubbleWriter!" );
      }

      streamWriter.toStream( output );
   }

   output << "\n</DataArray>\n";
   output << "</Points>\n";

   if ( mgr.write2D_ )
   {
      VTKMeshWriter::writeConnectivityP2TrianglesWithBubble( mgr.vtkDataFormat_, output, storage, level );
   }
   else
   {
      WALBERLA_ABORT( "Missing 3D implementation inside VTKP2PlusBubbleWriter!" );
   }

   output << "<PointData>\n";

   // write all scalar P2PlusBubbleFunctions of supported value type
   for ( const auto& function : mgr.feFunctionRegistry_.getP2PlusBubbleFunctions().getFunctions< double >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP2PlusBubbleFunctions().getFunctions< float >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP2PlusBubbleFunctions().getFunctions< int32_t >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP2PlusBubbleFunctions().getFunctions< int64_t >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }

   // write all P2PlusBubbleVectorFunctions of supported value type
   for ( const auto& function : mgr.feFunctionRegistry_.getP2PlusBubbleVectorFunctions().getFunctions< double >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP2PlusBubbleVectorFunctions().getFunctions< float >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP2PlusBubbleVectorFunctions().getFunctions< int32_t >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP2PlusBubbleVectorFunctions().getFunctions< int64_t >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }

   output << "</PointData>\n";

   vtk::writePieceFooter( output );
}

template < typename value_t >
void VTKP2PlusBubbleWriter::writeScalarFunction( std::ostream&                              output,
                                                 const P2PlusBubbleFunction< value_t >&     function,
                                                 const std::shared_ptr< PrimitiveStorage >& storage,
                                                 const uint_t&                              level,
                                                 bool                                       write2D,
                                                 vtk::DataFormat                            vtkDataFormat )

{
   WALBERLA_ASSERT_EQUAL( storage, function.getStorage() );

   vtk::openDataElement( output, typeToString< value_t >(), function.getFunctionName(), 1, vtkDataFormat );

   VTKStreamWriter< value_t >   streamWriter( vtkDataFormat );
   const P2Function< value_t >& p2Function = reinterpret_cast< const P2Function< value_t >& >( function );
   VTKP2Writer::writeP2FunctionData( write2D, streamWriter, p2Function, storage, level );

   if ( write2D )
   {
      for ( const auto& itFaces : storage->getFaces() )
      {
         const Face&        face   = *itFaces.second;
         const PrimitiveID& faceID = itFaces.first;

         value_t* vertexData = face.getData( function.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         value_t* edgeData   = face.getData( function.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

         auto       bubbleData = function.getVolumeDoFFunction().dofMemory( faceID, level );
         const auto memLayout  = function.getVolumeDoFFunction().memoryLayout();

         for ( const auto& faceType : facedof::allFaceTypes )
         {
            for ( const auto& microFace : facedof::macroface::Iterator( level, faceType, 0 ) )
            {
               // obtain data indices of dofs associated with vertices and edges of micro-face
               std::array< uint_t, 3 > vertexDoFIndices;
               vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, faceType, level, vertexDoFIndices );

               std::array< uint_t, 3 > edgeDoFIndices;
               edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( microFace, faceType, level, edgeDoFIndices );

               // assemble "correction term"
               value_t vertexDoFsum = value_t( 0 );
               vertexDoFsum += vertexData[vertexDoFIndices[0]];
               vertexDoFsum += vertexData[vertexDoFIndices[1]];
               vertexDoFsum += vertexData[vertexDoFIndices[2]];

               value_t edgeDoFsum = value_t( 0 );
               edgeDoFsum += edgeData[edgeDoFIndices[0]];
               edgeDoFsum += edgeData[edgeDoFIndices[1]];
               edgeDoFsum += edgeData[edgeDoFIndices[2]];

               value_t correction = ( value_t( 4 ) * edgeDoFsum - vertexDoFsum ) / value_t( 9 );

               // as the other six shape functions are not zero at the center, the bubble dof
               // represents an excess, for VTK output we need to add to this excess the weighted
               // sum of the other shape functions
               streamWriter << bubbleData[volumedofspace::indexing::index(
                                   microFace.x(), microFace.y(), faceType, 0, 1, level, memLayout )] +
                                   correction;
            }
         }
      }
   }
   else
   {
      WALBERLA_ABORT( "Missing 3D implementation inside VTKP2PlusBubbleWriter!" );
   }

   streamWriter.toStream( output );

   output << "\n</DataArray>\n";
}

template < typename value_t >
void VTKP2PlusBubbleWriter::writeVectorFunction( std::ostream&                                output,
                                                 const P2PlusBubbleVectorFunction< value_t >& vFunction,
                                                 const std::shared_ptr< PrimitiveStorage >&   storage,
                                                 const uint_t&                                level,
                                                 bool                                         write2D,
                                                 vtk::DataFormat                              vtkDataFormat )

{
   WALBERLA_ASSERT_EQUAL( storage, vFunction.getStorage() );

   uint_t dim = vFunction.getDimension();

   vtk::openDataElement( output, typeToString< value_t >(), vFunction.getFunctionName(), dim, vtkDataFormat );

   VTKStreamWriter< value_t >         streamWriter( vtkDataFormat );
   const P2VectorFunction< value_t >& p2VecFunction = reinterpret_cast< const P2VectorFunction< value_t >& >( vFunction );
   VTKP2Writer::writeP2VectorFunctionData( write2D, streamWriter, p2VecFunction, storage, level );

   if ( write2D )
   {
      for ( const auto& itFaces : storage->getFaces() )
      {
         const Face&        face   = *itFaces.second;
         const PrimitiveID& faceID = itFaces.first;

         for ( const auto& faceType : facedof::allFaceTypes )
         {
            for ( const auto& microFace : facedof::macroface::Iterator( level, faceType, 0 ) )
            {
               for ( uint_t idx = 0; idx < dim; ++idx )
               {
                  value_t* vertexData =
                      face.getData( vFunction[idx].getVertexDoFFunction().getFaceDataID() )->getPointer( level );
                  value_t* edgeData = face.getData( vFunction[idx].getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

                  auto       bubbleData = vFunction[idx].getVolumeDoFFunction().dofMemory( faceID, level );
                  const auto memLayout  = vFunction[idx].getVolumeDoFFunction().memoryLayout();

                  // obtain data indices of dofs associated with vertices and edges of micro-face
                  std::array< uint_t, 3 > vertexDoFIndices;
                  vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, faceType, level, vertexDoFIndices );

                  std::array< uint_t, 3 > edgeDoFIndices;
                  edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( microFace, faceType, level, edgeDoFIndices );

                  // assemble "correction term"
                  value_t vertexDoFsum = value_t( 0 );
                  vertexDoFsum += vertexData[vertexDoFIndices[0]];
                  vertexDoFsum += vertexData[vertexDoFIndices[1]];
                  vertexDoFsum += vertexData[vertexDoFIndices[2]];

                  value_t edgeDoFsum = value_t( 0 );
                  edgeDoFsum += edgeData[edgeDoFIndices[0]];
                  edgeDoFsum += edgeData[edgeDoFIndices[1]];
                  edgeDoFsum += edgeData[edgeDoFIndices[2]];

                  value_t correction = ( value_t( 4 ) * edgeDoFsum - vertexDoFsum ) / value_t( 9 );

                  // as the other six shape functions are not zero at the center, the bubble dof
                  // represents an excess, for VTK output we need to add to this excess the weighted
                  // sum of the other shape functions
                  streamWriter << bubbleData[volumedofspace::indexing::index(
                                      microFace.x(), microFace.y(), faceType, 0, 1, level, memLayout )] +
                                      correction;
               }
            }
         }
      }
   }
   else
   {
      WALBERLA_ABORT( "Missing 3D implementation inside VTKP2PlusBubbleWriter!" );
   }

   streamWriter.toStream( output );

   output << "\n</DataArray>\n";
}

} // namespace hyteg
