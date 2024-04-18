/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
                     const auto vtkPoint = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[i] );
                     Point3D    xBlend;
                     face.getGeometryMap()->evalF( vtkPoint, xBlend );
                     dstStream << xBlend[0] << xBlend[1] << xBlend[2];
                  }
               }
            }
         }
         else
         {
            size_t  rowsize = levelinfo::num_microvertices_per_edge( level );
            Point3D x, x0, xBlend;

            x0 = face.getCoordinates()[0];

            Point3D d0 = ( face.getCoordinates()[1] - face.getCoordinates()[0] ) / ( real_c( rowsize ) - 1 );
            Point3D d2 = ( face.getCoordinates()[2] - face.getCoordinates()[0] ) / ( real_c( rowsize ) - 1 );

            size_t inner_rowsize = rowsize;

            for ( size_t i = 0; i < rowsize; ++i )
            {
               x = x0;
               x += real_c( i ) * d2;

               for ( size_t j = 0; j < inner_rowsize; ++j )
               {
                  face.getGeometryMap()->evalF( x, xBlend );
                  dstStream << xBlend[0] << xBlend[1] << xBlend[2];
                  x += d0;
               }

               --inner_rowsize;
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
                     const Point3D vtkPoint = vertexdof::macrocell::coordinateFromIndex( level, cell, vIdx );
                     Point3D       xBlend;
                     cell.getGeometryMap()->evalF( vtkPoint, xBlend );
                     dstStream << xBlend[0] << xBlend[1] << xBlend[2];
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
               const Point3D vtkPoint = vertexdof::macrocell::coordinateFromIndex( level, cell, idxIt );
               Point3D       xBlend;
               cell.getGeometryMap()->evalF( vtkPoint, xBlend );
               dstStream << xBlend[0] << xBlend[1] << xBlend[2];
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

         const Point3D faceBottomLeftCoords  = face.getCoordinates()[0];
         const Point3D faceBottomRightCoords = face.getCoordinates()[1];
         const Point3D faceTopLeftCoords     = face.getCoordinates()[2];

         const Point3D horizontalMicroEdgeOffset =
             ( ( faceBottomRightCoords - faceBottomLeftCoords ) / real_c( levelinfo::num_microedges_per_edge( level ) ) ) * 0.5;
         const Point3D verticalMicroEdgeOffset =
             ( ( faceTopLeftCoords - faceBottomLeftCoords ) / real_c( levelinfo::num_microedges_per_edge( level ) ) ) * 0.5;

         Point3D xBlend;

         switch ( dofType )
         {
         case vtk::DoFType::EDGE_X: {
            for ( const auto& itIdx : edgedof::macroface::Iterator( level, 0 ) )
            {
               const Point3D horizontalMicroEdgePosition =
                   faceBottomLeftCoords + ( real_c( itIdx.x() * 2 + 1 ) * horizontalMicroEdgeOffset +
                                            real_c( itIdx.y() * 2 ) * verticalMicroEdgeOffset );
               face.getGeometryMap()->evalF( horizontalMicroEdgePosition, xBlend );
               dstStream << xBlend[0] << xBlend[1] << xBlend[2];
            }
            break;
         }
         case vtk::DoFType::EDGE_Y: {
            for ( const auto& itIdx : edgedof::macroface::Iterator( level, 0 ) )
            {
               const Point3D verticalMicroEdgePosition =
                   faceBottomLeftCoords + ( real_c( itIdx.x() * 2 ) * horizontalMicroEdgeOffset +
                                            real_c( itIdx.y() * 2 + 1 ) * verticalMicroEdgeOffset );
               face.getGeometryMap()->evalF( verticalMicroEdgePosition, xBlend );
               dstStream << xBlend[0] << xBlend[1] << xBlend[2];
            }
            break;
         }
         case vtk::DoFType::EDGE_XY: {
            for ( const auto& itIdx : edgedof::macroface::Iterator( level, 0 ) )
            {
               const Point3D horizontalMicroEdgePosition =
                   faceBottomLeftCoords + ( real_c( itIdx.x() * 2 + 1 ) * horizontalMicroEdgeOffset +
                                            real_c( itIdx.y() * 2 ) * verticalMicroEdgeOffset );
               const Point3D diagonalMicroEdgePosition = horizontalMicroEdgePosition + verticalMicroEdgeOffset;
               face.getGeometryMap()->evalF( diagonalMicroEdgePosition, xBlend );
               dstStream << xBlend[0] << xBlend[1] << xBlend[2];
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
               dstStream << microEdgePosition[0] << microEdgePosition[1] << microEdgePosition[2];
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
   VTKOutput::VTKStreamWriter< CellType > streamWriterCells( vtkDataFormat );
   writeElementNodeAssociationP1Triangles( streamWriterCells, storage, faceWidth, discontinuous );
   streamWriterCells.toStream( output );

   output << "\n</DataArray>\n";

   using OffsetType = uint32_t;

   vtk::openDataElement( output, typeToString< OffsetType >(), "offsets", 0, vtkDataFormat );

   VTKOutput::VTKStreamWriter< OffsetType > streamWriterOffsets( vtkDataFormat );

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

   VTKOutput::VTKStreamWriter< CellTypeType > streamWriterTypes( vtkDataFormat );

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
   // TODO: this should? be consistent with the type in the StreamWriterCell, but that performs a cast anyhow.
   using CellType  = uint32_t;
   CellType offset = 0;

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
         CellType rowsize       = static_cast< CellType >( faceWidth ) - 1;
         CellType inner_rowsize = rowsize;

         for ( CellType i = 0; i < rowsize; ++i )
         {
            for ( CellType j = 0; j < inner_rowsize - 1; ++j )
            {
               dstStream << offset << offset + 1u << offset + inner_rowsize + 1u;
               dstStream << offset + 1 << offset + inner_rowsize + 2u << offset + inner_rowsize + 1u;
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
      VTKOutput::VTKStreamWriter< CellType > streamWriterCells( vtkDataFormat );
      writeElementNodeAssociationP2Triangles( streamWriterCells, storage, level );
      streamWriterCells.toStream( output );
   }

   output << "\n</DataArray>\n";

   using OffsetType = uint32_t;

   vtk::openDataElement( output, typeToString< OffsetType >(), "offsets", 0, vtkDataFormat );

   VTKOutput::VTKStreamWriter< OffsetType > streamWriterOffsets( vtkDataFormat );

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

   VTKOutput::VTKStreamWriter< CellTypeType > streamWriterTypes( vtkDataFormat );

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

   using CellType  = uint32_t;
   CellType offset = 0;

   for ( auto& it : storage->getFaces() )
   {
      WALBERLA_UNUSED( it );

      // we execute the loops on a (virtually) refined mesh, so indices will fit to the vertices written
      // with writePointsForMicroVertices for P2 on (level+1)
      CellType rowsize       = static_cast< CellType >( levelinfo::num_microvertices_per_edge( level + 1 ) ) - 1;
      CellType inner_rowsize = rowsize;
      CellType idx0{}, idx1{}, idx2{}, idx3{}, idx4{}, idx5{};

      VTK_QUADRATIC_TRIANGLE_LOG( "rowsize = " << rowsize );

      for ( int i = 0; i <= static_cast< int >( rowsize ) - 3; i += 2 )
      {
         for ( int j = 0; j <= static_cast< int >( inner_rowsize ) - 4; j += 2 )
         {
            // lower left triangle
            idx0 = offset;
            idx1 = offset + 2 * inner_rowsize + 1;
            idx2 = offset + 2;
            idx3 = offset + inner_rowsize + 1;
            idx4 = offset + inner_rowsize + 2;
            idx5 = offset + 1;

            dstStream << idx0 << idx1 << idx2 << idx3 << idx4 << idx5;
            VTK_QUADRATIC_TRIANGLE_LOG( "[" << i << " , " << j << "]:" << std::setw( 2 ) << idx0 << " " << std::setw( 2 ) << idx1
                                            << " " << std::setw( 2 ) << idx2 << " " << std::setw( 2 ) << idx3 << " "
                                            << std::setw( 2 ) << idx4 << " " << std::setw( 2 ) << idx5 );

            // upper right triangle
            idx0 = offset + 2 * inner_rowsize + 1;
            idx1 = idx0 + 2;
            idx2 = offset + 2;
            idx3 = idx0 + 1;
            idx4 = offset + inner_rowsize + 3;
            idx5 = idx4 - 1;

            dstStream << idx0 << idx1 << idx2 << idx3 << idx4 << idx5;
            VTK_QUADRATIC_TRIANGLE_LOG( "[" << i << " , " << j << "]:" << std::setw( 2 ) << idx0 << " " << std::setw( 2 ) << idx1
                                            << " " << std::setw( 2 ) << idx2 << " " << std::setw( 2 ) << idx3 << " "
                                            << std::setw( 2 ) << idx4 << " " << std::setw( 2 ) << idx5 );

            // triangles live on refinement level "level"
            offset += 2;
         }

         // lower left triangle on top of column
         idx0 = offset;
         idx1 = offset + 2 * inner_rowsize + 1;
         idx2 = offset + 2;
         idx3 = offset + inner_rowsize + 1;
         idx4 = offset + inner_rowsize + 2;
         idx5 = offset + 1;

         dstStream << idx0 << idx1 << idx2 << idx3 << idx4 << idx5;
         VTK_QUADRATIC_TRIANGLE_LOG( "[" << i << " , *]:" << std::setw( 2 ) << idx0 << " " << std::setw( 2 ) << idx1 << " "
                                         << std::setw( 2 ) << idx2 << " " << std::setw( 2 ) << idx3 << " " << std::setw( 2 )
                                         << idx4 << " " << std::setw( 2 ) << idx5 );

         // skip on column w.r.t. fine level
         offset += 2 + inner_rowsize + 1;
         inner_rowsize -= 2;
      }

      // lower left triangle on the very "right"
      idx0 = offset;
      idx1 = offset + 2 * inner_rowsize + 1;
      idx2 = offset + 2;
      idx3 = offset + inner_rowsize + 1;
      idx4 = offset + inner_rowsize + 2;
      idx5 = offset + 1;

      dstStream << idx0 << idx1 << idx2 << idx3 << idx4 << idx5;
      VTK_QUADRATIC_TRIANGLE_LOG( "[* , *]:" << std::setw( 2 ) << idx0 << " " << std::setw( 2 ) << idx1 << " " << std::setw( 2 )
                                             << idx2 << " " << std::setw( 2 ) << idx3 << " " << std::setw( 2 ) << idx4 << " "
                                             << std::setw( 2 ) << idx5 );

      // prepare offset for next cell
      offset += 6;
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

   VTKOutput::VTKStreamWriter< CellIdx_T > streamWriterCells( vtkDataFormat );
   writeElementNodeAssociationP1Tetrahedrons( streamWriterCells, storage, width, discontinuous );
   streamWriterCells.toStream( output );

   output << "\n</DataArray>\n";

   using OffsetType = uint32_t;
   vtk::openDataElement( output, typeToString< OffsetType >(), "offsets", 0, vtkDataFormat );

   VTKOutput::VTKStreamWriter< OffsetType > streamWriterOffsets( vtkDataFormat );

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

   VTKOutput::VTKStreamWriter< CellTypeType > streamWriterTypes( vtkDataFormat );

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
   using CellIdx_T               = uint32_t;
   CellIdx_T offset              = 0;

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
               offset += 4;
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
      VTKOutput::VTKStreamWriter< CellIdx_T > streamWriterCells( vtkDataFormat );
      writeElementNodeAssociationP2Tetrahedrons( streamWriterCells, storage, level );
      streamWriterCells.toStream( output );
   }

   output << "\n</DataArray>\n";

   // this is the number of vtkQuadraticTetra cells in the mesh
   const uint_t numberOfCells = levelinfo::num_microcells_per_cell( level );

   using OffsetType           = uint32_t;
   vtk::openDataElement( output, typeToString< OffsetType >(), "offsets", 0, vtkDataFormat );

   VTKOutput::VTKStreamWriter< OffsetType > streamWriterOffsets( vtkDataFormat );

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

   VTKOutput::VTKStreamWriter< CellTypeType > streamWriterTypes( vtkDataFormat );

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

   // calculates the position of the point in the VTK list of points from a logical vertex index
   auto calcVTKPointArrayPosition = [level]( const indexing::Index& vertexIndex ) -> uint_t {
      const uint_t width{ levelinfo::num_microvertices_per_edge( level + 1 ) };
      const uint_t zOffset = levelinfo::num_microvertices_per_cell_from_width( width ) -
                             levelinfo::num_microvertices_per_cell_from_width( width - uint_c( vertexIndex.z() ) );
      const uint_t yOffset =
          levelinfo::num_microvertices_per_face_from_width( width - uint_c( vertexIndex.z() ) ) -
          levelinfo::num_microvertices_per_face_from_width( width - uint_c( vertexIndex.z() ) - uint_c( vertexIndex.y() ) );
      const uint_t xOffset = uint_c( vertexIndex.x() );
      return xOffset + yOffset + zOffset;
   };

   auto writeConnectivityForAllCellsOfType = [&dstStream, &edgeByVertexIndices, calcVTKPointArrayPosition](
                                                 const uint_t length, const uint_t baseIdx, celldof::CellType cellType ) {
      for ( const auto& it : indexing::CellIterator( length ) )
      {
         const auto spanningVertexIndices = celldof::macrocell::getMicroVerticesFromMicroCell( it, cellType );

         // write indices for the tetrahedrons vertices
         for ( const auto& spanningVertexIndex : spanningVertexIndices )
         {
            dstStream << baseIdx + calcVTKPointArrayPosition( 2 * spanningVertexIndex );
         }

         // write indices for the tetrahedrons edge midpoints
         for ( uint_t k = 0; k < 6; ++k )
         {
            celldof::Index vtxOne  = spanningVertexIndices[edgeByVertexIndices[k][0]];
            celldof::Index vtxTwo  = spanningVertexIndices[edgeByVertexIndices[k][1]];
            celldof::Index edgeIdx = ( vtxOne + vtxTwo );
            dstStream << baseIdx + calcVTKPointArrayPosition( edgeIdx );
         }
      }
   };

   // This is the number of (virtual) vertices in the mesh, i.e. the number of DoF locations
   const uint_t numberOfVertices =
       levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microvertices_per_edge( level + 1 ) );

   const uint_t width{ levelinfo::num_microvertices_per_edge( level ) };

   for ( uint_t macroCellIdx = 0; macroCellIdx < storage->getNumberOfLocalCells(); macroCellIdx++ )
   {
      uint_t baseIdx{ macroCellIdx * numberOfVertices };
      writeConnectivityForAllCellsOfType( width - 1, baseIdx, celldof::CellType::WHITE_UP );
      writeConnectivityForAllCellsOfType( width - 2, baseIdx, celldof::CellType::BLUE_UP );
      writeConnectivityForAllCellsOfType( width - 2, baseIdx, celldof::CellType::GREEN_UP );
      writeConnectivityForAllCellsOfType( width - 3, baseIdx, celldof::CellType::WHITE_DOWN );
      writeConnectivityForAllCellsOfType( width - 2, baseIdx, celldof::CellType::BLUE_DOWN );
      writeConnectivityForAllCellsOfType( width - 2, baseIdx, celldof::CellType::GREEN_DOWN );
   }
}

// ========================
//  explicit instantiation
// ========================
template void VTKMeshWriter::writePointsForMicroVertices( bool                                       write2D,
                                                          VTKOutput::VTKStreamWriter< real_t >&      dstStream,
                                                          const std::shared_ptr< PrimitiveStorage >& storage,
                                                          uint_t                                     level,
                                                          bool                                       discontinuous );

template void VTKMeshWriter::writePointsForMicroEdges( bool                                       write2D,
                                                       VTKOutput::VTKStreamWriter< real_t >&      dstStream,
                                                       const std::shared_ptr< PrimitiveStorage >& storage,
                                                       uint_t                                     level,
                                                       const vtk::DoFType&                        dofType );

#ifdef HYTEG_BUILD_WITH_ADIOS2
template void VTKMeshWriter::writePointsForMicroVertices( bool                                       write2D,
                                                          AdiosWriter::StreamAccessBuffer< real_t >& dstStream,
                                                          const std::shared_ptr< PrimitiveStorage >& storage,
                                                          uint_t                                     level,
                                                          bool                                       discontinuous );

template void VTKMeshWriter::writeElementNodeAssociationP1Triangles( AdiosWriter::StreamAccessBuffer< uint64_t, 4 >& dstStream,
                                                                     const std::shared_ptr< PrimitiveStorage >&      storage,
                                                                     uint_t                                          faceWidth,
                                                                     bool discontinuous );

template void VTKMeshWriter::writeElementNodeAssociationP1Tetrahedrons( AdiosWriter::StreamAccessBuffer< uint64_t, 5 >& dstStream,
                                                                        const std::shared_ptr< PrimitiveStorage >&      storage,
                                                                        uint_t                                          width,
                                                                        bool discontinuous );

template void VTKMeshWriter::writeElementNodeAssociationP2Triangles( AdiosWriter::StreamAccessBuffer< uint64_t, 7 >& dstStream_t,
                                                                     const std::shared_ptr< PrimitiveStorage >&      storage,
                                                                     uint_t                                          level );

template void
    VTKMeshWriter::writeElementNodeAssociationP2Tetrahedrons( AdiosWriter::StreamAccessBuffer< uint64_t, 11 >& dstStream_t,
                                                              const std::shared_ptr< PrimitiveStorage >&      storage,
                                                              uint_t                                          level );

#endif

} // namespace hyteg
