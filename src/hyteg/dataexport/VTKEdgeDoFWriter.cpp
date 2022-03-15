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
#include "hyteg/dataexport/VTKEdgeDoFWriter.hpp"

#include "core/DataTypes.h"

#include "hyteg/dataexport/VTKHelpers.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"

// from walberla
#include "vtk/UtilityFunctions.h"

namespace hyteg {

using walberla::vtk::typeToString;

void VTKEdgeDoFWriter::write( const VTKOutput& mgr, std::ostream& output, uint_t level, const vtk::DoFType& dofType )
{
   WALBERLA_ASSERT( dofType == vtk::DoFType::EDGE_X || dofType == vtk::DoFType::EDGE_Y || dofType == vtk::DoFType::EDGE_Z ||
                    dofType == vtk::DoFType::EDGE_XY || dofType == vtk::DoFType::EDGE_XZ || dofType == vtk::DoFType::EDGE_YZ ||
                    dofType == vtk::DoFType::EDGE_XYZ );

   if ( mgr.edgeDoFFunctions_.size() == 0 )
   {
      return;
   }

   auto storage = mgr.storage_;

   const uint_t numberOfPoints2D = storage->getNumberOfLocalFaces() * levelinfo::num_microedges_per_face( level ) / 3;
   const uint_t faceWidth        = levelinfo::num_microedges_per_edge( level );
   const uint_t numberOfCells2D  = storage->getNumberOfLocalFaces() * ( ( ( ( faceWidth - 1 ) * faceWidth ) / 2 ) +
                                                                       ( ( ( faceWidth - 2 ) * ( faceWidth - 1 ) ) / 2 ) );

   const uint_t numberOfPoints3D = storage->getNumberOfLocalCells() * [dofType, level]() {
      if ( dofType == vtk::DoFType::EDGE_XYZ )
         return levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) - 1 );
      else
         return levelinfo::num_microvertices_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) );
   }();
   const uint_t numberOfCells3D = storage->getNumberOfLocalCells() * [dofType, level]() {
      if ( dofType == vtk::DoFType::EDGE_XYZ )
         return levelinfo::num_microcells_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) - 1 );
      else
         return levelinfo::num_microcells_per_cell_from_width( levelinfo::num_microedges_per_edge( level ) );
   }();

   if ( mgr.write2D_ )
   {
      vtk::writePieceHeader( output, numberOfPoints2D, numberOfCells2D );
   }
   else
   {
      vtk::writePieceHeader( output, numberOfPoints3D, numberOfCells3D );
   }
   output << "<Points>\n";
   vtk::openDataElement( output, typeToString< real_t >(), "", 3, mgr.vtkDataFormat_ );

   VTKMeshWriter::writePointsForMicroEdges( mgr, output, storage, level, dofType );

   output << "\n</DataArray>\n";
   output << "</Points>\n";

   output << "<PointData>\n";

   for ( const auto& function : mgr.edgeDoFFunctions_.getFunctions< double >() )
   {
      writeScalarFunction( mgr, output, function, storage, level, dofType );
   }
   for ( const auto& function : mgr.edgeDoFFunctions_.getFunctions< int32_t >() )
   {
      writeScalarFunction( mgr, output, function, storage, level, dofType );
   }
   for ( const auto& function : mgr.edgeDoFFunctions_.getFunctions< int64_t >() )
   {
      writeScalarFunction( mgr, output, function, storage, level, dofType );
   }

   output << "</PointData>\n";

   if ( mgr.write2D_ )
   {
      VTKMeshWriter::writeCells2D( mgr, output, storage, levelinfo::num_microedges_per_edge( level ) );
   }
   else
   {
      if ( dofType == vtk::DoFType::EDGE_XYZ )
         VTKMeshWriter::writeCells3D( mgr, output, storage, levelinfo::num_microedges_per_edge( level ) - 1 );
      else
         VTKMeshWriter::writeCells3D( mgr, output, storage, levelinfo::num_microedges_per_edge( level ) );
   }

   vtk::writePieceFooter( output );
}

template < typename value_t >
void VTKEdgeDoFWriter::writeScalarFunction( const VTKOutput&                           mgr,
                                            std::ostream&                              output,
                                            const EdgeDoFFunction< value_t >&          function,
                                            const std::shared_ptr< PrimitiveStorage >& storage,
                                            uint_t                                     level,
                                            const vtk::DoFType&                        dofType )
{
   WALBERLA_ASSERT_EQUAL( storage, function.getStorage() );

   vtk::openDataElement( output, typeToString< value_t >(), function.getFunctionName(), 1, mgr.vtkDataFormat_ );

   WALBERLA_ASSERT( dofType == vtk::DoFType::EDGE_X || dofType == vtk::DoFType::EDGE_Y || dofType == vtk::DoFType::EDGE_Z ||
                    dofType == vtk::DoFType::EDGE_XY || dofType == vtk::DoFType::EDGE_XZ || dofType == vtk::DoFType::EDGE_YZ ||
                    dofType == vtk::DoFType::EDGE_XYZ );

   VTKOutput::VTKStreamWriter< value_t > streamWriter( mgr.vtkDataFormat_ );

   if ( mgr.write2D_ )
   {
      for ( const auto& it : storage->getFaces() )
      {
         const Face& face = *it.second;

         switch ( dofType )
         {
         case vtk::DoFType::EDGE_X: {
            for ( const auto& itIdx : edgedof::macroface::Iterator( level ) )
            {
               streamWriter << face.getData( function.getFaceDataID() )
                                   ->getPointer( level )[edgedof::macroface::horizontalIndex( level, uint_c( itIdx.col() ), uint_c( itIdx.row() ) )];
            }
            break;
         }
         case vtk::DoFType::EDGE_Y: {
            for ( const auto& itIdx : edgedof::macroface::Iterator( level ) )
            {
               streamWriter << face.getData( function.getFaceDataID() )
                                   ->getPointer( level )[edgedof::macroface::verticalIndex( level, uint_c( itIdx.col() ), uint_c( itIdx.row() ) )];
            }
            break;
         }
         case vtk::DoFType::EDGE_XY: {
            for ( const auto& itIdx : edgedof::macroface::Iterator( level ) )
            {
               streamWriter << face.getData( function.getFaceDataID() )
                                   ->getPointer( level )[edgedof::macroface::diagonalIndex( level, uint_c( itIdx.col() ), uint_c( itIdx.row() ) )];
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
      for ( const auto& it : storage->getCells() )
      {
         const Cell& cell     = *it.second;
         const auto  cellData = cell.getData( function.getCellDataID() )->getPointer( level );

         if ( dofType == vtk::DoFType::EDGE_XYZ )
         {
            for ( const auto& itIdx : edgedof::macrocell::IteratorXYZ( level ) )
            {
               streamWriter << cellData[edgedof::macrocell::xyzIndex( level, uint_c( itIdx.x() ), uint_c( itIdx.y() ), uint_c( itIdx.z() ) )];
            }
         }
         else
         {
            for ( const auto& itIdx : edgedof::macrocell::Iterator( level ) )
            {
               uint_t idx;
               switch ( dofType )
               {
               case vtk::DoFType::EDGE_X:
                  idx = edgedof::macrocell::xIndex( level, uint_c(itIdx.x() ), uint_c( itIdx.y() ), uint_c( itIdx.z() ) );
                  break;
               case vtk::DoFType::EDGE_Y:
                  idx = edgedof::macrocell::yIndex( level, uint_c( itIdx.x() ), uint_c( itIdx.y() ), uint_c( itIdx.z() ) );
                  break;
               case vtk::DoFType::EDGE_Z:
                  idx = edgedof::macrocell::zIndex( level, uint_c( itIdx.x() ), uint_c( itIdx.y() ), uint_c( itIdx.z() ) );
                  break;
               case vtk::DoFType::EDGE_XY:
                  idx = edgedof::macrocell::xyIndex( level, uint_c( itIdx.x() ), uint_c( itIdx.y() ), uint_c( itIdx.z() ) );
                  break;
               case vtk::DoFType::EDGE_XZ:
                  idx = edgedof::macrocell::xzIndex( level, uint_c( itIdx.x() ), uint_c( itIdx.y() ), uint_c( itIdx.z() ) );
                  break;
               case vtk::DoFType::EDGE_YZ:
                  idx = edgedof::macrocell::yzIndex( level, uint_c( itIdx.x() ), uint_c( itIdx.y() ), uint_c( itIdx.z() ) );
                  break;
               default:
                  WALBERLA_ABORT( "[VTK] Invalid DoFType" );
                  break;
               }
               streamWriter << cellData[idx];
            }
         }
      }
   }

   streamWriter.toStream( output );

   output << "\n</DataArray>\n";
}

} // namespace hyteg
