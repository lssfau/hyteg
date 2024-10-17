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
#include "core/DataTypes.h"

#include "hyteg/dataexport/VTKOutput/VTKHelpers.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/mesh/micro/MicroMesh.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"

#ifdef HYTEG_BUILD_WITH_ADIOS2
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#endif

// from walberla
#include "vtk/UtilityFunctions.h"

namespace hyteg {

using walberla::vtk::typeToString;

template < typename dstStream_t >
void VTKMeshWriter::writePointsForMicroVertices( bool                                       write2D,
                                                 dstStream_t&                               dstStream,
                                                 const std::shared_ptr< PrimitiveStorage >& storage,
                                                 uint_t                                     level,
                                                 bool                                       discontinuous )
{
   if ( write2D )
   {
      for ( const auto& it : storage->getFaces() )
      {
         Face& face = *it.second;

         if ( discontinuous )
         {
            for ( auto faceType : facedof::allFaceTypes )
            {
               for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
               {
                  const std::array< indexing::Index, 3 > vertexIndices =
                      facedof::macroface::getMicroVerticesFromMicroFace( idxIt, faceType );
                  for ( uint_t i = 0; i < 3; i++ )
                  {
                     const auto vtkPoint = micromesh::microVertexPosition( storage, face.getID(), level, vertexIndices[i] );
                     dstStream << vtkPoint[0] << vtkPoint[1] << vtkPoint[2];
                  }
               }
            }
         }
         else
         {
            for ( const auto& idxIt : vertexdof::macroface::Iterator( level, 0 ) )
            {
               const Point3D vtkPoint = micromesh::microVertexPosition( storage, face.getID(), level, idxIt );
               dstStream << vtkPoint[0] << vtkPoint[1] << vtkPoint[2];
            }
         }
      }
   }
   else
   {
      if ( discontinuous )
      {
         for ( const auto& it : storage->getCells() )
         {
            const Cell& cell = *it.second;

            for ( auto cellType : celldof::allCellTypes )
            {
               for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
               {
                  auto vertexIndices = celldof::macrocell::getMicroVerticesFromMicroCell( idxIt, cellType );
                  for ( auto vIdx : vertexIndices )
                  {
                     const auto vtkPoint = micromesh::microVertexPosition( storage, cell.getID(), level, vIdx );
                     dstStream << vtkPoint[0] << vtkPoint[1] << vtkPoint[2];
                  }
               }
            }
         }
      }
      else
      {
         for ( const auto& it : storage->getCells() )
         {
            const Cell& cell = *it.second;

            for ( const auto& idxIt : vertexdof::macrocell::Iterator( level, 0 ) )
            {
               const auto vtkPoint = micromesh::microVertexPosition( storage, cell.getID(), level, idxIt );
               dstStream << vtkPoint[0] << vtkPoint[1] << vtkPoint[2];
            }
         }
      }
   }
}

template < typename dstStream_t >
void VTKMeshWriter::writePointsForMicroEdges( bool                                       write2D,
                                              dstStream_t&                               dstStream,
                                              const std::shared_ptr< PrimitiveStorage >& storage,
                                              uint_t                                     level,
                                              const vtk::DoFType&                        dofType )
{
   if ( write2D )
   {
      WALBERLA_ASSERT( dofType == vtk::DoFType::EDGE_X || dofType == vtk::DoFType::EDGE_Y || dofType == vtk::DoFType::EDGE_XY );

      for ( const auto& it : storage->getFaces() )
      {
         Face& face = *it.second;

         switch ( dofType )
         {
         case vtk::DoFType::EDGE_X: {
            for ( const auto& itIdx : edgedof::macroface::Iterator( level, 0 ) )
            {
               const Point3D horizontalMicroEdgePosition =
                   micromesh::microEdgeCenterPosition( storage, face.getID(), level, itIdx, edgedof::EdgeDoFOrientation::X );
               dstStream << horizontalMicroEdgePosition[0] << horizontalMicroEdgePosition[1] << horizontalMicroEdgePosition[2];
            }
            break;
         }
         case vtk::DoFType::EDGE_Y: {
            for ( const auto& itIdx : edgedof::macroface::Iterator( level, 0 ) )
            {
               const Point3D verticalMicroEdgePosition =
                   micromesh::microEdgeCenterPosition( storage, face.getID(), level, itIdx, edgedof::EdgeDoFOrientation::Y );
               dstStream << verticalMicroEdgePosition[0] << verticalMicroEdgePosition[1] << verticalMicroEdgePosition[2];
            }
            break;
         }
         case vtk::DoFType::EDGE_XY: {
            for ( const auto& itIdx : edgedof::macroface::Iterator( level, 0 ) )
            {
               const Point3D diagonalMicroEdgePosition =
                   micromesh::microEdgeCenterPosition( storage, face.getID(), level, itIdx, edgedof::EdgeDoFOrientation::XY );
               dstStream << diagonalMicroEdgePosition[0] << diagonalMicroEdgePosition[1] << diagonalMicroEdgePosition[2];
            }
            break;
         }
         default:
            WALBERLA_ABORT( "Bad DoF type in VTK output for edge DoFs" );
            break;
         }
      }
   }
   else
   {
      WALBERLA_ASSERT( dofType == vtk::DoFType::EDGE_X || dofType == vtk::DoFType::EDGE_Y || dofType == vtk::DoFType::EDGE_Z ||
                       dofType == vtk::DoFType::EDGE_XY || dofType == vtk::DoFType::EDGE_XZ || dofType == vtk::DoFType::EDGE_YZ ||
                       dofType == vtk::DoFType::EDGE_XYZ );

      for ( const auto& it : storage->getCells() )
      {
         Cell&   cell = *it.second;
         Point3D microEdgePosition;

         if ( dofType == vtk::DoFType::EDGE_XYZ )
         {
            for ( const auto& itIdx : edgedof::macrocell::IteratorXYZ( level, 0 ) )
            {
               microEdgePosition =
                   micromesh::microEdgeCenterPosition( storage, cell.getID(), level, itIdx, edgedof::EdgeDoFOrientation::XYZ );
               dstStream << microEdgePosition[0] << microEdgePosition[1] << microEdgePosition[2];
            }
         }
         else
         {
            for ( const auto& itIdx : edgedof::macrocell::Iterator( level, 0 ) )
            {
               switch ( dofType )
               {
               case vtk::DoFType::EDGE_X:
                  microEdgePosition =
                      micromesh::microEdgeCenterPosition( storage, cell.getID(), level, itIdx, edgedof::EdgeDoFOrientation::X );
                  break;
               case vtk::DoFType::EDGE_Y:
                  microEdgePosition =
                      micromesh::microEdgeCenterPosition( storage, cell.getID(), level, itIdx, edgedof::EdgeDoFOrientation::Y );
                  break;
               case vtk::DoFType::EDGE_Z:
                  microEdgePosition =
                      micromesh::microEdgeCenterPosition( storage, cell.getID(), level, itIdx, edgedof::EdgeDoFOrientation::Z );
                  break;
               case vtk::DoFType::EDGE_XY:
                  microEdgePosition =
                      micromesh::microEdgeCenterPosition( storage, cell.getID(), level, itIdx, edgedof::EdgeDoFOrientation::XY );
                  break;
               case vtk::DoFType::EDGE_XZ:
                  microEdgePosition =
                      micromesh::microEdgeCenterPosition( storage, cell.getID(), level, itIdx, edgedof::EdgeDoFOrientation::XZ );
                  break;
               case vtk::DoFType::EDGE_YZ:
                  microEdgePosition =
                      micromesh::microEdgeCenterPosition( storage, cell.getID(), level, itIdx, edgedof::EdgeDoFOrientation::YZ );

                  break;
               default:
                  WALBERLA_ABORT( "[VTK] Invalid vtk::DoFType" );
                  break;
               }
               dstStream << microEdgePosition[0] << microEdgePosition[1] << microEdgePosition[2];
            }
         }
      }
   }
}

void VTKMeshWriter::writeCells2D( vtk::DataFormat                            vtkDataFormat,
                                  std::ostream&                              output,
                                  const std::shared_ptr< PrimitiveStorage >& storage,
                                  uint_t                                     faceWidth,
                                  bool                                       discontinuous )
{
   using CellType  = uint32_t;
   CellType offset = 0;

   output << "<Cells>\n";
   vtk::openDataElement( output, typeToString< CellType >(), "connectivity", 0, vtkDataFormat );

   const uint_t numberOfCells = levelinfo::num_microfaces_per_face_from_width( faceWidth );

   // connectivity
   VTKStreamWriter< CellType > streamWriterCells( vtkDataFormat );
   writeElementNodeAssociationP1Triangles( streamWriterCells, storage, faceWidth, discontinuous );
   streamWriterCells.toStream( output );

   output << "\n</DataArray>\n";

   using OffsetType = uint32_t;

   vtk::openDataElement( output, typeToString< OffsetType >(), "offsets", 0, vtkDataFormat );

   VTKStreamWriter< OffsetType > streamWriterOffsets( vtkDataFormat );

   // offsets
   offset = 3;
   for ( auto& it : storage->getFaces() )
   {
      WALBERLA_UNUSED( it );

      for ( OffsetType i = 0; i < numberOfCells; ++i )
      {
         streamWriterOffsets << offset;
         offset += 3;
      }
   }

   streamWriterOffsets.toStream( output );

   output << "\n</DataArray>\n";

   using CellTypeType = uint16_t;

   vtk::openDataElement( output, typeToString< CellTypeType >(), "types", 0, vtkDataFormat );

   VTKStreamWriter< CellTypeType > streamWriterTypes( vtkDataFormat );

   // cell types
   for ( auto& it : storage->getFaces() )
   {
      WALBERLA_UNUSED( it );
      for ( size_t i = 0; i < numberOfCells; ++i )
      {
         streamWriterTypes << 5;
      }
   }

   streamWriterTypes.toStream( output );

   output << "\n</DataArray>\n";
   output << "</Cells>\n";
}

template < typename dstStream_t >
void VTKMeshWriter::writeElementNodeAssociationP1Triangles( dstStream_t&                               dstStream,
                                                            const std::shared_ptr< PrimitiveStorage >& storage,
                                                            uint_t                                     faceWidth,
                                                            bool                                       discontinuous )
{
   uint_t offset = 0u;

   for ( auto& it : storage->getFaces() )
   {
      WALBERLA_UNUSED( it );

      if ( discontinuous )
      {
         for ( auto faceType : facedof::allFaceTypes )
         {
            const uint_t numMicroFacesAtBoundary = faceType == facedof::FaceType::GRAY ? faceWidth - 1 : faceWidth - 2;
            for ( const auto& idxIt : indexing::FaceIterator( numMicroFacesAtBoundary ) )
            {
               dstStream << offset << offset + 1u << offset + 2u;
               offset += 3u;
               WALBERLA_UNUSED( idxIt );
            }
         }
      }
      else
      {
         uint_t rowsize       = static_cast< uint_t >( faceWidth ) - 1u;
         uint_t inner_rowsize = rowsize;

         for ( uint_t i = 0; i < rowsize; ++i )
         {
            for ( uint_t j = 0; j < inner_rowsize - 1; ++j )
            {
               dstStream << offset << offset + 1u << offset + inner_rowsize + 1u;
               dstStream << offset + 1u << offset + inner_rowsize + 2u << offset + inner_rowsize + 1u;
               ++offset;
            }

            dstStream << offset << offset + 1u << offset + inner_rowsize + 1u;

            offset += 2u;
            --inner_rowsize;
         }

         // prepare offset for next cell
         // only necessary for continuous elements
         ++offset;
      }
   }
}

void VTKMeshWriter::writeConnectivityP2Triangles( vtk::DataFormat                            vtkDataFormat,
                                                  std::ostream&                              output,
                                                  const std::shared_ptr< PrimitiveStorage >& storage,
                                                  uint_t                                     level,
                                                  bool                                       discontinuous )
{
   using CellType = uint32_t;

   output << "<Cells>\n";
   vtk::openDataElement( output, typeToString< CellType >(), "connectivity", 0, vtkDataFormat );

   const uint_t numberOfCells = levelinfo::num_microfaces_per_face( level );

   if ( discontinuous )
   {
      WALBERLA_ABORT( "writeConnectivityP2Triangles does not support discontinous == true, yet!" );
   }

   else
   {
      VTKStreamWriter< CellType > streamWriterCells( vtkDataFormat );
      writeElementNodeAssociationP2Triangles( streamWriterCells, storage, level );
      streamWriterCells.toStream( output );
   }

   output << "\n</DataArray>\n";

   using OffsetType = uint32_t;

   vtk::openDataElement( output, typeToString< OffsetType >(), "offsets", 0, vtkDataFormat );

   VTKStreamWriter< OffsetType > streamWriterOffsets( vtkDataFormat );

   // offsets
   CellType offset = 6;
   for ( auto& it : storage->getFaces() )
   {
      WALBERLA_UNUSED( it );

      for ( OffsetType i = 0; i < numberOfCells; ++i )
      {
         streamWriterOffsets << offset;
         offset += 6;
      }
   }

   streamWriterOffsets.toStream( output );

   output << "\n</DataArray>\n";

   using CellTypeType = uint16_t;

   vtk::openDataElement( output, typeToString< CellTypeType >(), "types", 0, vtkDataFormat );

   VTKStreamWriter< CellTypeType > streamWriterTypes( vtkDataFormat );

   // cell types
   const unsigned char vtkQuadraticTriangleID = 22;

   for ( auto& it : storage->getFaces() )
   {
      WALBERLA_UNUSED( it );
      for ( size_t i = 0; i < numberOfCells; ++i )
      {
         streamWriterTypes << vtkQuadraticTriangleID;
      }
   }

   streamWriterTypes.toStream( output );

   output << "\n</DataArray>\n";
   output << "</Cells>\n";
}

template < typename dstStream_t >
void VTKMeshWriter::writeElementNodeAssociationP2Triangles( dstStream_t&                               dstStream,
                                                            const std::shared_ptr< PrimitiveStorage >& storage,
                                                            uint_t                                     level )
{
// #define DEBUG_VTK_QUADRATIC_TRIANGLE
#ifdef DEBUG_VTK_QUADRATIC_TRIANGLE
#define VTK_QUADRATIC_TRIANGLE_LOG( msg ) WALBERLA_LOG_INFO_ON_ROOT( msg );
#else
#define VTK_QUADRATIC_TRIANGLE_LOG( msg )
#endif

   uint_t macroFaceCount = 0;

   for ( auto& it : storage->getFaces() )
   {
      WALBERLA_UNUSED( it );

      VTK_QUADRATIC_TRIANGLE_LOG( "rowsize = " << rowsize );

      std::map< vtk::DoFType, uint_t > offsets;

      offsets[vtk::DoFType::VERTEX] = macroFaceCount * levelinfo::num_microvertices_per_face( level );
      offsets[vtk::DoFType::EDGE_X] = storage->getNumberOfLocalFaces() * levelinfo::num_microvertices_per_face( level ) +
                                      macroFaceCount * ( levelinfo::num_microedges_per_face( level ) / 3 );
      offsets[vtk::DoFType::EDGE_XY] = storage->getNumberOfLocalFaces() * ( levelinfo::num_microvertices_per_face( level ) +
                                                                            levelinfo::num_microedges_per_face( level ) / 3 ) +
                                       macroFaceCount * ( levelinfo::num_microedges_per_face( level ) / 3 );
      offsets[vtk::DoFType::EDGE_Y] =
          storage->getNumberOfLocalFaces() *
              ( levelinfo::num_microvertices_per_face( level ) + 2 * ( levelinfo::num_microedges_per_face( level ) / 3 ) ) +
          macroFaceCount * ( levelinfo::num_microedges_per_face( level ) / 3 );

      for ( auto faceType : facedof::allFaceTypes )
      {
         for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
         {
            if ( faceType == facedof::FaceType::GRAY )
            {
               uint_t idx0 = offsets[vtk::DoFType::VERTEX] + vertexdof::macroface::index( level, idxIt.x(), idxIt.y() );
               uint_t idx1 = offsets[vtk::DoFType::VERTEX] + vertexdof::macroface::index( level, idxIt.x() + 1, idxIt.y() );
               uint_t idx2 = offsets[vtk::DoFType::VERTEX] + vertexdof::macroface::index( level, idxIt.x(), idxIt.y() + 1 );

               // Looks odd, but we need the index of the x-edge for all types due to the ordering of the points.
               // The index offsets are included in the offsets map.

               uint_t idx3 = offsets[vtk::DoFType::EDGE_X] +
                             edgedof::macroface::index( level, idxIt.x(), idxIt.y(), edgedof::EdgeDoFOrientation::X );

               uint_t idx4 = offsets[vtk::DoFType::EDGE_XY] +
                             edgedof::macroface::index( level, idxIt.x(), idxIt.y(), edgedof::EdgeDoFOrientation::X );

               uint_t idx5 = offsets[vtk::DoFType::EDGE_Y] +
                             edgedof::macroface::index( level, idxIt.x(), idxIt.y(), edgedof::EdgeDoFOrientation::X );

               dstStream << idx0 << idx1 << idx2 << idx3 << idx4 << idx5;
            }

            if ( faceType == facedof::FaceType::BLUE )
            {
               uint_t idx0 = offsets[vtk::DoFType::VERTEX] + vertexdof::macroface::index( level, idxIt.x() + 1, idxIt.y() );
               uint_t idx1 = offsets[vtk::DoFType::VERTEX] + vertexdof::macroface::index( level, idxIt.x() + 1, idxIt.y() + 1 );
               uint_t idx2 = offsets[vtk::DoFType::VERTEX] + vertexdof::macroface::index( level, idxIt.x(), idxIt.y() + 1 );

               uint_t idx3 = offsets[vtk::DoFType::EDGE_Y] +
                             edgedof::macroface::index( level, idxIt.x() + 1, idxIt.y(), edgedof::EdgeDoFOrientation::X );

               uint_t idx4 = offsets[vtk::DoFType::EDGE_X] +
                             edgedof::macroface::index( level, idxIt.x(), idxIt.y() + 1, edgedof::EdgeDoFOrientation::X );

               uint_t idx5 = offsets[vtk::DoFType::EDGE_XY] +
                             edgedof::macroface::index( level, idxIt.x(), idxIt.y(), edgedof::EdgeDoFOrientation::X );

               dstStream << idx0 << idx1 << idx2 << idx3 << idx4 << idx5;
            }
         }
      }

      macroFaceCount++;
   }
}

void VTKMeshWriter::writeCells3D( vtk::DataFormat                            vtkDataFormat,
                                  std::ostream&                              output,
                                  const std::shared_ptr< PrimitiveStorage >& storage,
                                  uint_t                                     width,
                                  bool                                       discontinuous )
{
   using CellIdx_T            = int32_t;
   CellIdx_T    offset        = 0;
   const uint_t numberOfCells = levelinfo::num_microcells_per_cell_from_width( width );

   output << "<Cells>\n";
   vtk::openDataElement( output, typeToString< CellIdx_T >(), "connectivity", 0, vtkDataFormat );

   VTKStreamWriter< CellIdx_T > streamWriterCells( vtkDataFormat );
   writeElementNodeAssociationP1Tetrahedrons( streamWriterCells, storage, width, discontinuous );
   streamWriterCells.toStream( output );

   output << "\n</DataArray>\n";

   using OffsetType = uint32_t;
   vtk::openDataElement( output, typeToString< OffsetType >(), "offsets", 0, vtkDataFormat );

   VTKStreamWriter< OffsetType > streamWriterOffsets( vtkDataFormat );

   // offsets
   offset = 4;
   for ( const auto& it : storage->getCells() )
   {
      WALBERLA_UNUSED( it );

      for ( size_t i = 0; i < numberOfCells; ++i )
      {
         streamWriterOffsets << offset;
         offset += 4;
      }
   }

   streamWriterOffsets.toStream( output );

   output << "\n</DataArray>\n";

   using CellTypeType = uint16_t;
   vtk::openDataElement( output, typeToString< CellTypeType >(), "types", 0, vtkDataFormat );

   VTKStreamWriter< CellTypeType > streamWriterTypes( vtkDataFormat );

   // cell types
   for ( const auto& it : storage->getCells() )
   {
      WALBERLA_UNUSED( it );
      for ( size_t i = 0; i < numberOfCells; ++i )
      {
         streamWriterTypes << 10;
      }
   }

   streamWriterTypes.toStream( output );

   output << "\n</DataArray>\n";
   output << "</Cells>\n";
}

template < typename dstStream_t >
void VTKMeshWriter::writeElementNodeAssociationP1Tetrahedrons( dstStream_t&                               dstStream,
                                                               const std::shared_ptr< PrimitiveStorage >& storage,
                                                               uint_t                                     width,
                                                               bool                                       discontinuous )
{
   // calculates the position of the point in the VTK list of points from a logical vertex index
   auto calcVTKPointArrayPosition = [width]( const indexing::Index& vertexIndex ) -> uint_t {
      const uint_t zOffset = levelinfo::num_microvertices_per_cell_from_width( width ) -
                             levelinfo::num_microvertices_per_cell_from_width( width - uint_c( vertexIndex.z() ) );
      const uint_t yOffset =
          levelinfo::num_microvertices_per_face_from_width( width - uint_c( vertexIndex.z() ) ) -
          levelinfo::num_microvertices_per_face_from_width( width - uint_c( vertexIndex.z() ) - uint_c( vertexIndex.y() ) );
      const uint_t xOffset = uint_c( vertexIndex.x() );
      return xOffset + yOffset + zOffset;
   };

   const uint_t numberOfVertices = levelinfo::num_microvertices_per_cell_from_width( width );
   uint_t       offset           = 0u;

   for ( uint_t macroCellIdx = 0; macroCellIdx < storage->getNumberOfLocalCells(); macroCellIdx++ )
   {
      if ( discontinuous )
      {
         for ( auto cellType : celldof::allCellTypes )
         {
            celldof::macrocell::numCellsPerRowByTypeFromWidth( width, cellType );
            const uint_t numMicroCellsAtBoundary = celldof::macrocell::numCellsPerRowByTypeFromWidth( width, cellType );
            for ( const auto& idxIt : indexing::CellIterator( numMicroCellsAtBoundary ) )
            {
               dstStream << offset << offset + 1u << offset + 2u << offset + 3u;
               offset += 4u;
               WALBERLA_UNUSED( idxIt );
            }
         }
      }
      else
      {
         for ( const auto& it : indexing::CellIterator( width - 1 ) )
         {
            const auto spanningVertexIndices =
                celldof::macrocell::getMicroVerticesFromMicroCell( it, celldof::CellType::WHITE_UP );

            for ( const auto& spanningVertexIndex : spanningVertexIndices )
            {
               dstStream << macroCellIdx * numberOfVertices + calcVTKPointArrayPosition( spanningVertexIndex );
            }
         }

         for ( const auto& it : indexing::CellIterator( width - 2 ) )
         {
            const auto spanningVertexIndices =
                celldof::macrocell::getMicroVerticesFromMicroCell( it, celldof::CellType::BLUE_UP );

            for ( const auto& spanningVertexIndex : spanningVertexIndices )
            {
               dstStream << macroCellIdx * numberOfVertices + calcVTKPointArrayPosition( spanningVertexIndex );
            }
         }

         for ( const auto& it : indexing::CellIterator( width - 2 ) )
         {
            const auto spanningVertexIndices =
                celldof::macrocell::getMicroVerticesFromMicroCell( it, celldof::CellType::GREEN_UP );

            for ( const auto& spanningVertexIndex : spanningVertexIndices )
            {
               dstStream << macroCellIdx * numberOfVertices + calcVTKPointArrayPosition( spanningVertexIndex );
            }
         }

         for ( const auto& it : indexing::CellIterator( width - 3 ) )
         {
            const auto spanningVertexIndices =
                celldof::macrocell::getMicroVerticesFromMicroCell( it, celldof::CellType::WHITE_DOWN );

            for ( const auto& spanningVertexIndex : spanningVertexIndices )
            {
               dstStream << macroCellIdx * numberOfVertices + calcVTKPointArrayPosition( spanningVertexIndex );
            }
         }

         for ( const auto& it : indexing::CellIterator( width - 2 ) )
         {
            const auto spanningVertexIndices =
                celldof::macrocell::getMicroVerticesFromMicroCell( it, celldof::CellType::BLUE_DOWN );

            for ( const auto& spanningVertexIndex : spanningVertexIndices )
            {
               dstStream << macroCellIdx * numberOfVertices + calcVTKPointArrayPosition( spanningVertexIndex );
            }
         }

         for ( const auto& it : indexing::CellIterator( width - 2 ) )
         {
            const auto spanningVertexIndices =
                celldof::macrocell::getMicroVerticesFromMicroCell( it, celldof::CellType::GREEN_DOWN );

            for ( const auto& spanningVertexIndex : spanningVertexIndices )
            {
               dstStream << macroCellIdx * numberOfVertices + calcVTKPointArrayPosition( spanningVertexIndex );
            }
         }
      }
   }
}

void VTKMeshWriter::writeConnectivityP2Tetrahedrons( vtk::DataFormat                            vtkDataFormat,
                                                     std::ostream&                              output,
                                                     const std::shared_ptr< PrimitiveStorage >& storage,
                                                     uint_t                                     level,
                                                     bool                                       discontinuous )
{
   using CellIdx_T = int32_t;

   if ( discontinuous )
   {
      WALBERLA_ABORT( "writeConnectivityP2Tetrahedrons does not support discontinous == true, yet!" );
   }

   output << "<Cells>\n";
   vtk::openDataElement( output, typeToString< CellIdx_T >(), "connectivity", 0, vtkDataFormat );

   {
      VTKStreamWriter< CellIdx_T > streamWriterCells( vtkDataFormat );
      writeElementNodeAssociationP2Tetrahedrons( streamWriterCells, storage, level );
      streamWriterCells.toStream( output );
   }

   output << "\n</DataArray>\n";

   // this is the number of vtkQuadraticTetra cells in the mesh
   const uint_t numberOfCells = levelinfo::num_microcells_per_cell( level );

   using OffsetType = uint32_t;
   vtk::openDataElement( output, typeToString< OffsetType >(), "offsets", 0, vtkDataFormat );

   VTKStreamWriter< OffsetType > streamWriterOffsets( vtkDataFormat );

   // offsets
   uint_t offset = 10;
   for ( const auto& it : storage->getCells() )
   {
      WALBERLA_UNUSED( it );

      for ( size_t i = 0; i < numberOfCells; ++i )
      {
         streamWriterOffsets << offset;
         offset += 10;
      }
   }

   streamWriterOffsets.toStream( output );

   output << "\n</DataArray>\n";

   using CellTypeType = uint16_t;
   vtk::openDataElement( output, typeToString< CellTypeType >(), "types", 0, vtkDataFormat );

   VTKStreamWriter< CellTypeType > streamWriterTypes( vtkDataFormat );

   // cell types
   const unsigned char vtkQuadraticTetraID = 24;
   for ( const auto& it : storage->getCells() )
   {
      WALBERLA_UNUSED( it );
      for ( size_t i = 0; i < numberOfCells; ++i )
      {
         streamWriterTypes << vtkQuadraticTetraID;
      }
   }

   streamWriterTypes.toStream( output );

   output << "\n</DataArray>\n";
   output << "</Cells>\n";
}

template < typename dstStream_t >
void VTKMeshWriter::writeElementNodeAssociationP2Tetrahedrons( dstStream_t&                               dstStream,
                                                               const std::shared_ptr< PrimitiveStorage >& storage,
                                                               uint_t                                     level )
{
   // From the vtk documentation:
   //
   // vtkQuadraticTetra is a concrete implementation of vtkNonLinearCell to represent a three-dimensional, 10-node,
   // isoparametric parabolic tetrahedron. The interpolation is the standard finite element, quadratic isoparametric
   // shape function. The cell includes a mid-edge node on each of the size edges of the tetrahedron. The ordering
   // of the ten points defining the cell is point ids (0-3,4-9) where ids 0-3 are the four tetra vertices;
   // and point ids 4-9 are the midedge nodes between (0,1), (1,2), (2,0), (0,3), (1,3), and (2,3).
   const std::array< std::array< uint_t, 2 >, 6 > edgeByVertexIndices{ 0, 1, 1, 2, 2, 0, 0, 3, 1, 3, 2, 3 };

   const auto macroCells = storage->getNumberOfLocalCells();

   // calculates the position of the point in the VTK list of points from a logical vertex index
   auto calcVTKPointArrayPosition = [level, macroCells]( const indexing::Index& vertexIndex, uint_t macroCellIndex ) -> uint_t {
#if 0
      const uint_t width{ levelinfo::num_microvertices_per_edge( level + 1 ) };
      const uint_t zOffset = levelinfo::num_microvertices_per_cell_from_width( width ) -
                             levelinfo::num_microvertices_per_cell_from_width( width - uint_c( vertexIndex.z() ) );
      const uint_t yOffset =
          levelinfo::num_microvertices_per_face_from_width( width - uint_c( vertexIndex.z() ) ) -
          levelinfo::num_microvertices_per_face_from_width( width - uint_c( vertexIndex.z() ) - uint_c( vertexIndex.y() ) );
      const uint_t xOffset = uint_c( vertexIndex.x() );
      return xOffset + yOffset + zOffset;
#endif
      const auto x = vertexIndex.x();
      const auto y = vertexIndex.y();
      const auto z = vertexIndex.z();

      const auto xeven = x % 2 == 0;
      const auto yeven = y % 2 == 0;
      const auto zeven = z % 2 == 0;

      const auto xodd = !xeven;
      const auto yodd = !yeven;
      const auto zodd = !zeven;

      const auto numVertices    = levelinfo::num_microvertices_per_cell( level );
      const auto numEdgesButXYZ = levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) );
      const auto numEdgesXYZ =
          levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) - 1 );

      std::map< vtk::DoFType, uint_t > offsets;
      offsets[vtk::DoFType::VERTEX]  = macroCellIndex * numVertices;
      offsets[vtk::DoFType::EDGE_X]  = macroCells * ( numVertices + 0 * numEdgesButXYZ ) + macroCellIndex * numEdgesButXYZ;
      offsets[vtk::DoFType::EDGE_Y]  = macroCells * ( numVertices + 1 * numEdgesButXYZ ) + macroCellIndex * numEdgesButXYZ;
      offsets[vtk::DoFType::EDGE_Z]  = macroCells * ( numVertices + 2 * numEdgesButXYZ ) + macroCellIndex * numEdgesButXYZ;
      offsets[vtk::DoFType::EDGE_XY] = macroCells * ( numVertices + 3 * numEdgesButXYZ ) + macroCellIndex * numEdgesButXYZ;
      offsets[vtk::DoFType::EDGE_XZ] = macroCells * ( numVertices + 4 * numEdgesButXYZ ) + macroCellIndex * numEdgesButXYZ;
      offsets[vtk::DoFType::EDGE_YZ] = macroCells * ( numVertices + 5 * numEdgesButXYZ ) + macroCellIndex * numEdgesButXYZ;

      // Since we have no better indexing function, I am subtracting the offset of the xyz indexing function here.
      offsets[vtk::DoFType::EDGE_XYZ] = macroCells * ( numVertices + 6 * numEdgesButXYZ ) + macroCellIndex * numEdgesXYZ -
                                        edgedof::macrocell::index( level, 0, 0, 0, edgedof::EdgeDoFOrientation::XYZ );

      if ( xeven && yeven && zeven )
      {
         // vertexodof
         const auto idx = offsets[vtk::DoFType::VERTEX] + vertexdof::macrocell::index( level, x / 2, y / 2, z / 2 );
         return idx;
      }
      else
      {
         // edgedof

         // Looks odd (pun not intended), but we need the index of the x-edge for all types due to the ordering of the points.
         // The index offsets are included in the offsets map.

         if ( xodd && yeven && zeven )
         {
            const auto idx = offsets[vtk::DoFType::EDGE_X] +
                             edgedof::macrocell::index( level, ( x - 1 ) / 2, y / 2, z / 2, edgedof::EdgeDoFOrientation::X );
            return idx;
         }
         else if ( xeven && yodd && zeven )
         {
            const auto idx = offsets[vtk::DoFType::EDGE_Y] +
                             edgedof::macrocell::index( level, x / 2, ( y - 1 ) / 2, z / 2, edgedof::EdgeDoFOrientation::X );
            return idx;
         }
         else if ( xeven && yeven && zodd )
         {
            const auto idx = offsets[vtk::DoFType::EDGE_Z] +
                             edgedof::macrocell::index( level, x / 2, y / 2, ( z - 1 ) / 2, edgedof::EdgeDoFOrientation::X );
            return idx;
         }
         else if ( xodd && yodd && zeven )
         {
            const auto idx =
                offsets[vtk::DoFType::EDGE_XY] +
                edgedof::macrocell::index( level, ( x - 1 ) / 2, ( y - 1 ) / 2, z / 2, edgedof::EdgeDoFOrientation::X );
            return idx;
         }
         else if ( xodd && yeven && zodd )
         {
            const auto idx =
                offsets[vtk::DoFType::EDGE_XZ] +
                edgedof::macrocell::index( level, ( x - 1 ) / 2, y / 2, ( z - 1 ) / 2, edgedof::EdgeDoFOrientation::X );
            return idx;
         }
         else if ( xeven && yodd && zodd )
         {
            const auto idx =
                offsets[vtk::DoFType::EDGE_YZ] +
                edgedof::macrocell::index( level, x / 2, ( y - 1 ) / 2, ( z - 1 ) / 2, edgedof::EdgeDoFOrientation::X );
            return idx;
         }
         else if ( xodd && yodd && zodd )
         {
            const auto idx =
                offsets[vtk::DoFType::EDGE_XYZ] +
                edgedof::macrocell::index( level, ( x - 1 ) / 2, ( y - 1 ) / 2, ( z - 1 ) / 2, edgedof::EdgeDoFOrientation::XYZ );
            return idx;
         }
         else
         {
            WALBERLA_ABORT( "This should not happen." )
         }
      }
   };

   auto writeConnectivityForAllCellsOfType = [&dstStream, &edgeByVertexIndices, calcVTKPointArrayPosition](
                                                 const uint_t length, const uint_t macroCellIdx, celldof::CellType cellType ) {
      for ( const auto& it : indexing::CellIterator( length ) )
      {
         const auto spanningVertexIndices = celldof::macrocell::getMicroVerticesFromMicroCell( it, cellType );

         // write indices for the tetrahedrons vertices
         for ( const auto& spanningVertexIndex : spanningVertexIndices )
         {
            dstStream << calcVTKPointArrayPosition( 2 * spanningVertexIndex, macroCellIdx );
         }

         // write indices for the tetrahedrons edge midpoints
         for ( uint_t k = 0; k < 6; ++k )
         {
            celldof::Index vtxOne  = spanningVertexIndices[edgeByVertexIndices[k][0]];
            celldof::Index vtxTwo  = spanningVertexIndices[edgeByVertexIndices[k][1]];
            celldof::Index edgeIdx = ( vtxOne + vtxTwo );
            dstStream << calcVTKPointArrayPosition( edgeIdx, macroCellIdx );
         }
      }
   };

   const uint_t width{ levelinfo::num_microvertices_per_edge( level ) };

   for ( uint_t macroCellIdx = 0; macroCellIdx < storage->getNumberOfLocalCells(); macroCellIdx++ )
   {
      writeConnectivityForAllCellsOfType( width - 1u, macroCellIdx, celldof::CellType::WHITE_UP );
      writeConnectivityForAllCellsOfType( width - 2u, macroCellIdx, celldof::CellType::BLUE_UP );
      writeConnectivityForAllCellsOfType( width - 2u, macroCellIdx, celldof::CellType::GREEN_UP );
      writeConnectivityForAllCellsOfType( width - 3u, macroCellIdx, celldof::CellType::WHITE_DOWN );
      writeConnectivityForAllCellsOfType( width - 2u, macroCellIdx, celldof::CellType::BLUE_DOWN );
      writeConnectivityForAllCellsOfType( width - 2u, macroCellIdx, celldof::CellType::GREEN_DOWN );
   }
}

// ========================
//  explicit instantiation
// ========================
template void VTKMeshWriter::writePointsForMicroVertices( bool                                       write2D,
                                                          VTKStreamWriter< real_t >&                 dstStream,
                                                          const std::shared_ptr< PrimitiveStorage >& storage,
                                                          uint_t                                     level,
                                                          bool                                       discontinuous );

template void VTKMeshWriter::writePointsForMicroEdges( bool                                       write2D,
                                                       VTKStreamWriter< real_t >&                 dstStream,
                                                       const std::shared_ptr< PrimitiveStorage >& storage,
                                                       uint_t                                     level,
                                                       const vtk::DoFType&                        dofType );

#ifdef HYTEG_BUILD_WITH_ADIOS2
template void VTKMeshWriter::writePointsForMicroVertices( bool                                       write2D,
                                                          AdiosWriter::StreamAccessBuffer< real_t >& dstStream,
                                                          const std::shared_ptr< PrimitiveStorage >& storage,
                                                          uint_t                                     level,
                                                          bool                                       discontinuous );

template void VTKMeshWriter::writeElementNodeAssociationP1Triangles(
    AdiosWriter::StreamAccessBuffer< ADIOS2_PARAVIEW_INT_TYPE, 4 >& dstStream,
    const std::shared_ptr< PrimitiveStorage >&                      storage,
    uint_t                                                          faceWidth,
    bool                                                            discontinuous );

template void VTKMeshWriter::writeElementNodeAssociationP1Tetrahedrons(
    AdiosWriter::StreamAccessBuffer< ADIOS2_PARAVIEW_INT_TYPE, 5 >& dstStream,
    const std::shared_ptr< PrimitiveStorage >&                      storage,
    uint_t                                                          width,
    bool                                                            discontinuous );

template void VTKMeshWriter::writeElementNodeAssociationP2Triangles(
    AdiosWriter::StreamAccessBuffer< ADIOS2_PARAVIEW_INT_TYPE, 7 >& dstStream_t,
    const std::shared_ptr< PrimitiveStorage >&                      storage,
    uint_t                                                          level );

template void VTKMeshWriter::writeElementNodeAssociationP2Tetrahedrons(
    AdiosWriter::StreamAccessBuffer< ADIOS2_PARAVIEW_INT_TYPE, 11 >& dstStream_t,
    const std::shared_ptr< PrimitiveStorage >&                       storage,
    uint_t                                                           level );

template void VTKMeshWriter::writePointsForMicroEdges( bool                                       write2D,
                                                       AdiosWriter::StreamAccessBuffer< real_t >& dstStream,
                                                       const std::shared_ptr< PrimitiveStorage >& storage,
                                                       uint_t                                     level,
                                                       const vtk::DoFType&                        dofType );
#endif

} // namespace hyteg
