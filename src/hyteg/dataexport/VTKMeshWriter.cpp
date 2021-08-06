/*
 * Copyright (c) 2017-2021 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/dataexport/VTKDGDoFWriter.hpp"
#include "hyteg/dataexport/VTKHelpers.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"

// from walberla
#include "vtk/UtilityFunctions.h"

namespace hyteg {

using walberla::vtk::typeToString;

void VTKMeshWriter::writePointsForMicroVertices( const VTKOutput&                           mgr,
                                                 std::ostream&                              output,
                                                 const std::shared_ptr< PrimitiveStorage >& storage,
                                                 uint_t                                     level )
{
   using ScalarType = real_t;
   VTKOutput::VTKStreamWriter< ScalarType > streamWriter( mgr.vtkDataFormat_ );

   if ( mgr.write2D_ )
   {
      for ( const auto& it : storage->getFaces() )
      {
         Face& face = *it.second;

         size_t  rowsize = levelinfo::num_microvertices_per_edge( level );
         Point3D x, x0, xBlend;

         x0 = face.coords[0];

         Point3D d0 = ( face.coords[1] - face.coords[0] ) / ( real_c( rowsize ) - 1 );
         Point3D d2 = ( face.coords[2] - face.coords[0] ) / ( real_c( rowsize ) - 1 );

         size_t inner_rowsize = rowsize;

         for ( size_t i = 0; i < rowsize; ++i )
         {
            x = x0;
            x += real_c( i ) * d2;

            for ( size_t j = 0; j < inner_rowsize; ++j )
            {
               face.getGeometryMap()->evalF( x, xBlend );
               streamWriter << xBlend[0] << xBlend[1] << xBlend[2];
               x += d0;
            }

            --inner_rowsize;
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
            const Point3D vtkPoint = vertexdof::macrocell::coordinateFromIndex( level, cell, idxIt );
            Point3D       xBlend;
            cell.getGeometryMap()->evalF( vtkPoint, xBlend );
            streamWriter << xBlend[0] << xBlend[1] << xBlend[2];
         }
      }
   }

   streamWriter.toStream( output );
}

void VTKMeshWriter::writePointsForMicroEdges( const VTKOutput&                           mgr,
                                              std::ostream&                              output,
                                              const std::shared_ptr< PrimitiveStorage >& storage,
                                              uint_t                                     level,
                                              const vtk::DoFType&                        dofType )
{
   using ScalarType = real_t;
   VTKOutput::VTKStreamWriter< ScalarType > streamWriter( mgr.vtkDataFormat_ );

   if ( mgr.write2D_ )
   {
      WALBERLA_ASSERT( dofType == vtk::DoFType::EDGE_X || dofType == vtk::DoFType::EDGE_Y || dofType == vtk::DoFType::EDGE_XY );

      for ( const auto& it : storage->getFaces() )
      {
         Face& face = *it.second;

         const Point3D faceBottomLeftCoords  = face.coords[0];
         const Point3D faceBottomRightCoords = face.coords[1];
         const Point3D faceTopLeftCoords     = face.coords[2];

         const Point3D horizontalMicroEdgeOffset =
             ( ( faceBottomRightCoords - faceBottomLeftCoords ) / real_c( levelinfo::num_microedges_per_edge( level ) ) ) * 0.5;
         const Point3D verticalMicroEdgeOffset =
             ( ( faceTopLeftCoords - faceBottomLeftCoords ) / real_c( levelinfo::num_microedges_per_edge( level ) ) ) * 0.5;

         Point3D xBlend;

         switch ( dofType )
         {
         case vtk::DoFType::EDGE_X:
         {
            for ( const auto& itIdx : edgedof::macroface::Iterator( level, 0 ) )
            {
               const Point3D horizontalMicroEdgePosition =
                   faceBottomLeftCoords + ( real_c( itIdx.col() * 2 + 1 ) * horizontalMicroEdgeOffset +
                                            real_c( itIdx.row() * 2 ) * verticalMicroEdgeOffset );
               face.getGeometryMap()->evalF( horizontalMicroEdgePosition, xBlend );
               streamWriter << xBlend[0] << xBlend[1] << xBlend[2];
            }
            break;
         }
         case vtk::DoFType::EDGE_Y:
         {
            for ( const auto& itIdx : edgedof::macroface::Iterator( level, 0 ) )
            {
               const Point3D verticalMicroEdgePosition =
                   faceBottomLeftCoords + ( real_c( itIdx.col() * 2 ) * horizontalMicroEdgeOffset +
                                            real_c( itIdx.row() * 2 + 1 ) * verticalMicroEdgeOffset );
               face.getGeometryMap()->evalF( verticalMicroEdgePosition, xBlend );
               streamWriter << xBlend[0] << xBlend[1] << xBlend[2];
            }
            break;
         }
         case vtk::DoFType::EDGE_XY:
         {
            for ( const auto& itIdx : edgedof::macroface::Iterator( level, 0 ) )
            {
               const Point3D horizontalMicroEdgePosition =
                   faceBottomLeftCoords + ( real_c( itIdx.col() * 2 + 1 ) * horizontalMicroEdgeOffset +
                                            real_c( itIdx.row() * 2 ) * verticalMicroEdgeOffset );
               const Point3D diagonalMicroEdgePosition = horizontalMicroEdgePosition + verticalMicroEdgeOffset;
               face.getGeometryMap()->evalF( diagonalMicroEdgePosition, xBlend );
               streamWriter << xBlend[0] << xBlend[1] << xBlend[2];
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
               microEdgePosition = vertexdof::macrocell::coordinateFromIndex( level, cell, itIdx ) +
                                   edgedof::macrocell::xShiftFromVertex( level, cell ) +
                                   edgedof::macrocell::yShiftFromVertex( level, cell ) +
                                   edgedof::macrocell::zShiftFromVertex( level, cell );
               streamWriter << microEdgePosition[0] << microEdgePosition[1] << microEdgePosition[2];
            }
         }
         else
         {
            for ( const auto& itIdx : edgedof::macrocell::Iterator( level, 0 ) )
            {
               microEdgePosition = vertexdof::macrocell::coordinateFromIndex( level, cell, itIdx );
               switch ( dofType )
               {
               case vtk::DoFType::EDGE_X:
                  microEdgePosition += edgedof::macrocell::xShiftFromVertex( level, cell );
                  break;
               case vtk::DoFType::EDGE_Y:
                  microEdgePosition += edgedof::macrocell::yShiftFromVertex( level, cell );
                  break;
               case vtk::DoFType::EDGE_Z:
                  microEdgePosition += edgedof::macrocell::zShiftFromVertex( level, cell );
                  break;
               case vtk::DoFType::EDGE_XY:
                  microEdgePosition +=
                      edgedof::macrocell::xShiftFromVertex( level, cell ) + edgedof::macrocell::yShiftFromVertex( level, cell );
                  break;
               case vtk::DoFType::EDGE_XZ:
                  microEdgePosition +=
                      edgedof::macrocell::xShiftFromVertex( level, cell ) + edgedof::macrocell::zShiftFromVertex( level, cell );
                  break;
               case vtk::DoFType::EDGE_YZ:
                  microEdgePosition +=
                      edgedof::macrocell::yShiftFromVertex( level, cell ) + edgedof::macrocell::zShiftFromVertex( level, cell );
                  break;
               default:
                  WALBERLA_ABORT( "[VTK] Invalid vtk::DoFType" );
                  break;
               }
               streamWriter << microEdgePosition[0] << microEdgePosition[1] << microEdgePosition[2];
            }
         }
      }
   }

   streamWriter.toStream( output );
}

void VTKMeshWriter::writeCells2D( const VTKOutput&                           mgr,
                                  std::ostream&                              output,
                                  const std::shared_ptr< PrimitiveStorage >& storage,
                                  uint_t                                     faceWidth )
{
   using CellType = uint32_t;

   output << "<Cells>\n";
   vtk::openDataElement( output, typeToString< CellType >(), "connectivity", 0, mgr.vtkDataFormat_ );

   VTKOutput::VTKStreamWriter< CellType > streamWriterCells( mgr.vtkDataFormat_ );

   const uint_t numberOfCells =
       uint_c( ( ( faceWidth - 1 ) * faceWidth ) / 2 ) + ( ( ( faceWidth - 2 ) * ( faceWidth - 1 ) ) / 2 );

   // connectivity
   CellType offset = 0;

   for ( auto& it : storage->getFaces() )
   {
      //TODO is it really unused?
      WALBERLA_UNUSED( it );
      CellType rowsize       = static_cast< CellType >( faceWidth ) - 1;
      CellType inner_rowsize = rowsize;

      for ( CellType i = 0; i < rowsize; ++i )
      {
         for ( CellType j = 0; j < inner_rowsize - 1; ++j )
         {
            streamWriterCells << offset << offset + 1 << offset + inner_rowsize + 1;
            streamWriterCells << offset + 1 << offset + inner_rowsize + 2 << offset + inner_rowsize + 1;
            ++offset;
         }

         streamWriterCells << offset << offset + 1 << offset + inner_rowsize + 1;

         offset += 2;
         --inner_rowsize;
      }

      ++offset;
   }

   streamWriterCells.toStream( output );

   output << "\n</DataArray>\n";

   using OffsetType = uint32_t;

   vtk::openDataElement( output, typeToString< OffsetType >(), "offsets", 0, mgr.vtkDataFormat_ );

   VTKOutput::VTKStreamWriter< OffsetType > streamWriterOffsets( mgr.vtkDataFormat_ );

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

   vtk::openDataElement( output, typeToString< CellTypeType >(), "types", 0, mgr.vtkDataFormat_ );

   VTKOutput::VTKStreamWriter< CellTypeType > streamWriterTypes( mgr.vtkDataFormat_ );

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

void VTKMeshWriter::writeCells3D( const VTKOutput&                           mgr,
                                  std::ostream&                              output,
                                  const std::shared_ptr< PrimitiveStorage >& storage,
                                  uint_t                                     width )
{
   using CellIdx_T = int32_t;

   output << "<Cells>\n";
   vtk::openDataElement( output, typeToString< CellIdx_T >(), "connectivity", 0, mgr.vtkDataFormat_ );

   VTKOutput::VTKStreamWriter< CellIdx_T > streamWriterCells( mgr.vtkDataFormat_ );

   // calculates the position of the point in the VTK list of points from a logical vertex index
   auto calcVTKPointArrayPosition = [width]( const indexing::Index& vertexIndex ) -> uint_t {
      const uint_t zOffset = levelinfo::num_microvertices_per_cell_from_width( width ) -
                             levelinfo::num_microvertices_per_cell_from_width( width - vertexIndex.z() );
      const uint_t yOffset = levelinfo::num_microvertices_per_face_from_width( width - vertexIndex.z() ) -
                             levelinfo::num_microvertices_per_face_from_width( width - vertexIndex.z() - vertexIndex.y() );
      const uint_t xOffset = vertexIndex.x();
      return xOffset + yOffset + zOffset;
   };

   const uint_t numberOfVertices = levelinfo::num_microvertices_per_cell_from_width( width );
   const uint_t numberOfCells    = levelinfo::num_microcells_per_cell_from_width( width );

   for ( uint_t macroCellIdx = 0; macroCellIdx < storage->getNumberOfLocalCells(); macroCellIdx++ )
   {
      for ( const auto& it : indexing::CellIterator( width - 1 ) )
      {
         const auto spanningVertexIndices = celldof::macrocell::getMicroVerticesFromMicroCell( it, celldof::CellType::WHITE_UP );

         for ( const auto& spanningVertexIndex : spanningVertexIndices )
         {
            streamWriterCells << macroCellIdx * numberOfVertices + calcVTKPointArrayPosition( spanningVertexIndex );
         }
      }

      for ( const auto& it : indexing::CellIterator( width - 2 ) )
      {
         const auto spanningVertexIndices = celldof::macrocell::getMicroVerticesFromMicroCell( it, celldof::CellType::BLUE_UP );

         for ( const auto& spanningVertexIndex : spanningVertexIndices )
         {
            streamWriterCells << macroCellIdx * numberOfVertices + calcVTKPointArrayPosition( spanningVertexIndex );
         }
      }

      for ( const auto& it : indexing::CellIterator( width - 2 ) )
      {
         const auto spanningVertexIndices = celldof::macrocell::getMicroVerticesFromMicroCell( it, celldof::CellType::GREEN_UP );

         for ( const auto& spanningVertexIndex : spanningVertexIndices )
         {
            streamWriterCells << macroCellIdx * numberOfVertices + calcVTKPointArrayPosition( spanningVertexIndex );
         }
      }

      for ( const auto& it : indexing::CellIterator( width - 3 ) )
      {
         const auto spanningVertexIndices =
             celldof::macrocell::getMicroVerticesFromMicroCell( it, celldof::CellType::WHITE_DOWN );

         for ( const auto& spanningVertexIndex : spanningVertexIndices )
         {
            streamWriterCells << macroCellIdx * numberOfVertices + calcVTKPointArrayPosition( spanningVertexIndex );
         }
      }

      for ( const auto& it : indexing::CellIterator( width - 2 ) )
      {
         const auto spanningVertexIndices = celldof::macrocell::getMicroVerticesFromMicroCell( it, celldof::CellType::BLUE_DOWN );

         for ( const auto& spanningVertexIndex : spanningVertexIndices )
         {
            streamWriterCells << macroCellIdx * numberOfVertices + calcVTKPointArrayPosition( spanningVertexIndex );
         }
      }

      for ( const auto& it : indexing::CellIterator( width - 2 ) )
      {
         const auto spanningVertexIndices =
             celldof::macrocell::getMicroVerticesFromMicroCell( it, celldof::CellType::GREEN_DOWN );

         for ( const auto& spanningVertexIndex : spanningVertexIndices )
         {
            streamWriterCells << macroCellIdx * numberOfVertices + calcVTKPointArrayPosition( spanningVertexIndex );
         }
      }
   }

   streamWriterCells.toStream( output );

   output << "\n</DataArray>\n";

   using OffsetType = uint32_t;
   vtk::openDataElement( output, typeToString< OffsetType >(), "offsets", 0, mgr.vtkDataFormat_ );

   VTKOutput::VTKStreamWriter< OffsetType > streamWriterOffsets( mgr.vtkDataFormat_ );

   // offsets
   uint_t offset = 4;
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
   vtk::openDataElement( output, typeToString< CellTypeType >(), "types", 0, mgr.vtkDataFormat_ );

   VTKOutput::VTKStreamWriter< CellTypeType > streamWriterTypes( mgr.vtkDataFormat_ );

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

} // namespace hyteg
