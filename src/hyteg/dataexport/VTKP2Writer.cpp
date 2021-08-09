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

#include "hyteg/dataexport/VTKP2Writer.hpp"

#include "core/DataTypes.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/edgedofspace/EdgeDoFIndexing.hpp"

#include "vtk/UtilityFunctions.h"

namespace hyteg {

using walberla::vtk::typeToString;

void VTKP2Writer::write( const VTKOutput& mgr, std::ostream& output, const uint_t& level )
{
   if ( mgr.getNumRegisteredFunctions( vtk::DoFType::P2 ) == 0 )
   {
      return;
   }

   // auto storage = p2Functions_[0].getStorage();
   auto storage = mgr.storage_;

   const uint_t numberOfPoints = mgr.write2D_ ?
                                     storage->getNumberOfLocalFaces() * levelinfo::num_microvertices_per_face( level + 1 ) :
                                     storage->getNumberOfLocalCells() * levelinfo::num_microvertices_per_cell( level + 1 );
   const uint_t numberOfCells = mgr.write2D_ ?
                                    storage->getNumberOfLocalFaces() * levelinfo::num_microfaces_per_face( level + 1 ) :
                                    storage->getNumberOfLocalCells() * levelinfo::num_microcells_per_cell( level + 1 );
   ;

   vtk::writePieceHeader( output, numberOfPoints, numberOfCells );

   output << "<Points>\n";
   vtk::openDataElement( output, typeToString< real_t >(), "", 3, mgr.vtkDataFormat_ );

   VTKMeshWriter::writePointsForMicroVertices( mgr, output, storage, level + 1 );

   output << "\n</DataArray>\n";
   output << "</Points>\n";

   if ( mgr.write2D_ )
   {
      VTKMeshWriter::writeCells2D( mgr, output, storage, levelinfo::num_microvertices_per_edge( level + 1 ) );
   }
   else
   {
      VTKMeshWriter::writeCells3D( mgr, output, storage, levelinfo::num_microvertices_per_edge( level + 1 ) );
   }

   output << "<PointData>\n";

   for ( const auto& function : mgr.p2Functions_ )
   {
      writeScalarFunction( function, output, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }

#ifdef P2_VEC_SCALAR_OUTPUT
   for ( const auto& function : mgr.p2VecFunctions_ )
   {
      for ( uint_t idx = 0; idx < ( write2D_ ? 2 : 3 ); ++idx )
      {
         writeScalarFunction( function[idx], output, level );
      }
   }
#else
   for ( const auto& function : mgr.p2VecFunctions_ )
   {
      writeVectorFunction( function, output, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
#endif

   output << "</PointData>\n";

   vtk::writePieceFooter( output );
}

void VTKP2Writer::writeScalarFunction( const P2Function< real_t >& function,
                                       std::ostream&               output,
                                       const uint_t&               level,
                                       bool                        write2D,
                                       vtk::DataFormat             vtkDataFormat )

{
   using ScalarType = real_t;
   auto storage     = function.getStorage();
   vtk::openDataElement( output, typeToString< ScalarType >(), function.getFunctionName(), 1, vtkDataFormat );

   VTKOutput::VTKStreamWriter< ScalarType > streamWriter( vtkDataFormat );

   if ( write2D )
   {
      for ( const auto& itFaces : storage->getFaces() )
      {
         const Face& face = *itFaces.second;

         for ( const auto& it : vertexdof::macroface::Iterator( level + 1, 0 ) )
         {
            if ( it.row() % 2 == 0 )
            {
               if ( it.col() % 2 == 0 )
               {
                  streamWriter << face.getData( function.getVertexDoFFunction().getFaceDataID() )
                                      ->getPointer( level )[vertexdof::macroface::indexFromVertex(
                                          level, it.col() / 2, it.row() / 2, stencilDirection::VERTEX_C )];
               }
               else
               {
                  streamWriter
                      << face.getData( function.getEdgeDoFFunction().getFaceDataID() )
                             ->getPointer(
                                 level )[edgedof::macroface::horizontalIndex( level, ( it.col() - 1 ) / 2, it.row() / 2 )];
               }
            }
            else
            {
               if ( it.col() % 2 == 0 )
               {
                  streamWriter << face.getData( function.getEdgeDoFFunction().getFaceDataID() )
                                      ->getPointer(
                                          level )[edgedof::macroface::verticalIndex( level, it.col() / 2, ( it.row() - 1 ) / 2 )];
               }
               else
               {
                  streamWriter
                      << face.getData( function.getEdgeDoFFunction().getFaceDataID() )
                             ->getPointer(
                                 level )[edgedof::macroface::diagonalIndex( level, ( it.col() - 1 ) / 2, ( it.row() - 1 ) / 2 )];
               }
            }
         }
      }
   }
   else
   {
      for ( const auto& itCells : storage->getCells() )
      {
         const Cell& cell       = *itCells.second;
         auto        vertexData = cell.getData( function.getVertexDoFFunction().getCellDataID() )->getPointer( level );
         auto        edgeData   = cell.getData( function.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

         for ( const auto& it : vertexdof::macrocell::Iterator( level + 1, 0 ) )
         {
            const auto   x   = it.x();
            const auto   y   = it.y();
            const auto   z   = it.z();
            const uint_t mod = ( z % 2 << 0 ) | ( y % 2 << 1 ) | ( x % 2 << 2 );

            switch ( mod )
            {
            case 0b000:
               streamWriter
                   << vertexData[vertexdof::macrocell::indexFromVertex( level, x / 2, y / 2, z / 2, stencilDirection::VERTEX_C )];
               break;
            case 0b100:
               streamWriter << edgeData[edgedof::macrocell::xIndex( level, ( x - 1 ) / 2, y / 2, z / 2 )];
               break;
            case 0b010:
               streamWriter << edgeData[edgedof::macrocell::yIndex( level, x / 2, ( y - 1 ) / 2, z / 2 )];
               break;
            case 0b001:
               streamWriter << edgeData[edgedof::macrocell::zIndex( level, x / 2, y / 2, ( z - 1 ) / 2 )];
               break;
            case 0b110:
               streamWriter << edgeData[edgedof::macrocell::xyIndex( level, ( x - 1 ) / 2, ( y - 1 ) / 2, z / 2 )];
               break;
            case 0b101:
               streamWriter << edgeData[edgedof::macrocell::xzIndex( level, ( x - 1 ) / 2, y / 2, ( z - 1 ) / 2 )];
               break;
            case 0b011:
               streamWriter << edgeData[edgedof::macrocell::yzIndex( level, x / 2, ( y - 1 ) / 2, ( z - 1 ) / 2 )];
               break;
            case 0b111:
               streamWriter << edgeData[edgedof::macrocell::xyzIndex( level, ( x - 1 ) / 2, ( y - 1 ) / 2, ( z - 1 ) / 2 )];
               break;
            }
         }
      }
   }

   streamWriter.toStream( output );

   output << "\n</DataArray>\n";
}

void VTKP2Writer::writeVectorFunction( const P2VectorFunction< real_t >& function,
                                       std::ostream&                     output,
                                       const uint_t&                     level,
                                       bool                              write2D,
                                       vtk::DataFormat                   vtkDataFormat )
{
   using ScalarType = real_t;

   auto storage = function.getStorage();
   vtk::openDataElement(
       output, typeToString< ScalarType >(), function.getFunctionName(), function.getDimension(), vtkDataFormat );

   VTKOutput::VTKStreamWriter< ScalarType > streamWriter( vtkDataFormat );

   if ( write2D )
   {
      for ( const auto& itFaces : storage->getFaces() )
      {
         const Face& face = *itFaces.second;

         for ( const auto& it : vertexdof::macroface::Iterator( level + 1, 0 ) )
         {
            if ( it.row() % 2 == 0 )
            {
               if ( it.col() % 2 == 0 )
               {
                  streamWriter << face.getData( function[0].getVertexDoFFunction().getFaceDataID() )
                                      ->getPointer( level )[vertexdof::macroface::indexFromVertex(
                                          level, it.col() / 2, it.row() / 2, stencilDirection::VERTEX_C )];
                  streamWriter << face.getData( function[1].getVertexDoFFunction().getFaceDataID() )
                                      ->getPointer( level )[vertexdof::macroface::indexFromVertex(
                                          level, it.col() / 2, it.row() / 2, stencilDirection::VERTEX_C )];
               }
               else
               {
                  streamWriter
                      << face.getData( function[0].getEdgeDoFFunction().getFaceDataID() )
                             ->getPointer(
                                 level )[edgedof::macroface::horizontalIndex( level, ( it.col() - 1 ) / 2, it.row() / 2 )];
                  streamWriter
                      << face.getData( function[1].getEdgeDoFFunction().getFaceDataID() )
                             ->getPointer(
                                 level )[edgedof::macroface::horizontalIndex( level, ( it.col() - 1 ) / 2, it.row() / 2 )];
               }
            }
            else
            {
               if ( it.col() % 2 == 0 )
               {
                  streamWriter << face.getData( function[0].getEdgeDoFFunction().getFaceDataID() )
                                      ->getPointer(
                                          level )[edgedof::macroface::verticalIndex( level, it.col() / 2, ( it.row() - 1 ) / 2 )];
                  streamWriter << face.getData( function[1].getEdgeDoFFunction().getFaceDataID() )
                                      ->getPointer(
                                          level )[edgedof::macroface::verticalIndex( level, it.col() / 2, ( it.row() - 1 ) / 2 )];
               }
               else
               {
                  streamWriter
                      << face.getData( function[0].getEdgeDoFFunction().getFaceDataID() )
                             ->getPointer(
                                 level )[edgedof::macroface::diagonalIndex( level, ( it.col() - 1 ) / 2, ( it.row() - 1 ) / 2 )];
                  streamWriter
                      << face.getData( function[1].getEdgeDoFFunction().getFaceDataID() )
                             ->getPointer(
                                 level )[edgedof::macroface::diagonalIndex( level, ( it.col() - 1 ) / 2, ( it.row() - 1 ) / 2 )];
               }
            }
         }
      }
   }
   else
   {
      for ( const auto& itCells : storage->getCells() )
      {
         const Cell& cell = *itCells.second;

         auto vertexData0 = cell.getData( function[0].getVertexDoFFunction().getCellDataID() )->getPointer( level );
         auto vertexData1 = cell.getData( function[1].getVertexDoFFunction().getCellDataID() )->getPointer( level );
         auto vertexData2 = cell.getData( function[2].getVertexDoFFunction().getCellDataID() )->getPointer( level );

         auto edgeData0 = cell.getData( function[0].getEdgeDoFFunction().getCellDataID() )->getPointer( level );
         auto edgeData1 = cell.getData( function[1].getEdgeDoFFunction().getCellDataID() )->getPointer( level );
         auto edgeData2 = cell.getData( function[2].getEdgeDoFFunction().getCellDataID() )->getPointer( level );

         for ( const auto& it : vertexdof::macrocell::Iterator( level + 1, 0 ) )
         {
            const auto   x   = it.x();
            const auto   y   = it.y();
            const auto   z   = it.z();
            const uint_t mod = ( z % 2 << 0 ) | ( y % 2 << 1 ) | ( x % 2 << 2 );

            switch ( mod )
            {
            case 0b000:
               streamWriter << vertexData0[vertexdof::macrocell::indexFromVertex(
                   level, x / 2, y / 2, z / 2, stencilDirection::VERTEX_C )];
               streamWriter << vertexData1[vertexdof::macrocell::indexFromVertex(
                   level, x / 2, y / 2, z / 2, stencilDirection::VERTEX_C )];
               streamWriter << vertexData2[vertexdof::macrocell::indexFromVertex(
                   level, x / 2, y / 2, z / 2, stencilDirection::VERTEX_C )];
               break;
            case 0b100:
               streamWriter << edgeData0[edgedof::macrocell::xIndex( level, ( x - 1 ) / 2, y / 2, z / 2 )];
               streamWriter << edgeData1[edgedof::macrocell::xIndex( level, ( x - 1 ) / 2, y / 2, z / 2 )];
               streamWriter << edgeData2[edgedof::macrocell::xIndex( level, ( x - 1 ) / 2, y / 2, z / 2 )];
               break;
            case 0b010:
               streamWriter << edgeData0[edgedof::macrocell::yIndex( level, x / 2, ( y - 1 ) / 2, z / 2 )];
               streamWriter << edgeData1[edgedof::macrocell::yIndex( level, x / 2, ( y - 1 ) / 2, z / 2 )];
               streamWriter << edgeData2[edgedof::macrocell::yIndex( level, x / 2, ( y - 1 ) / 2, z / 2 )];
               break;
            case 0b001:
               streamWriter << edgeData0[edgedof::macrocell::zIndex( level, x / 2, y / 2, ( z - 1 ) / 2 )];
               streamWriter << edgeData1[edgedof::macrocell::zIndex( level, x / 2, y / 2, ( z - 1 ) / 2 )];
               streamWriter << edgeData2[edgedof::macrocell::zIndex( level, x / 2, y / 2, ( z - 1 ) / 2 )];
               break;
            case 0b110:
               streamWriter << edgeData0[edgedof::macrocell::xyIndex( level, ( x - 1 ) / 2, ( y - 1 ) / 2, z / 2 )];
               streamWriter << edgeData1[edgedof::macrocell::xyIndex( level, ( x - 1 ) / 2, ( y - 1 ) / 2, z / 2 )];
               streamWriter << edgeData2[edgedof::macrocell::xyIndex( level, ( x - 1 ) / 2, ( y - 1 ) / 2, z / 2 )];
               break;
            case 0b101:
               streamWriter << edgeData0[edgedof::macrocell::xzIndex( level, ( x - 1 ) / 2, y / 2, ( z - 1 ) / 2 )];
               streamWriter << edgeData1[edgedof::macrocell::xzIndex( level, ( x - 1 ) / 2, y / 2, ( z - 1 ) / 2 )];
               streamWriter << edgeData2[edgedof::macrocell::xzIndex( level, ( x - 1 ) / 2, y / 2, ( z - 1 ) / 2 )];
               break;
            case 0b011:
               streamWriter << edgeData0[edgedof::macrocell::yzIndex( level, x / 2, ( y - 1 ) / 2, ( z - 1 ) / 2 )];
               streamWriter << edgeData1[edgedof::macrocell::yzIndex( level, x / 2, ( y - 1 ) / 2, ( z - 1 ) / 2 )];
               streamWriter << edgeData2[edgedof::macrocell::yzIndex( level, x / 2, ( y - 1 ) / 2, ( z - 1 ) / 2 )];
               break;
            case 0b111:
               streamWriter << edgeData0[edgedof::macrocell::xyzIndex( level, ( x - 1 ) / 2, ( y - 1 ) / 2, ( z - 1 ) / 2 )];
               streamWriter << edgeData1[edgedof::macrocell::xyzIndex( level, ( x - 1 ) / 2, ( y - 1 ) / 2, ( z - 1 ) / 2 )];
               streamWriter << edgeData2[edgedof::macrocell::xyzIndex( level, ( x - 1 ) / 2, ( y - 1 ) / 2, ( z - 1 ) / 2 )];
               break;
            }
         }
      }
   }

   streamWriter.toStream( output );

   output << "\n</DataArray>\n";
}

} // namespace hyteg
