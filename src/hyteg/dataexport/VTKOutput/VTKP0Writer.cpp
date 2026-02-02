/*
 * Copyright (c) 2017-2026 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "hyteg/dataexport/VTKOutput/VTKP1Writer.hpp"

// from walberla
#include "vtk/UtilityFunctions.h"

namespace hyteg {

using walberla::vtk::typeToString;

void VTKP0Writer::write( const VTKOutput& mgr, std::ostream& output, const uint_t& level )
{
   if ( mgr.getNumRegisteredFunctions( vtk::DoFType::P0 ) == 0 )
   {
      return;
   }

   auto storage = mgr.storage_;

   const uint_t numberOfPoints2D = storage->getNumberOfLocalFaces() * levelinfo::num_microvertices_per_face( level );
   const uint_t numberOfCells2D  = storage->getNumberOfLocalFaces() * levelinfo::num_microfaces_per_face( level );

   const uint_t numberOfPoints3D = storage->getNumberOfLocalCells() * levelinfo::num_microvertices_per_cell( level );
   const uint_t numberOfCells3D  = storage->getNumberOfLocalCells() * levelinfo::num_microcells_per_cell( level );

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

   {
      VTKStreamWriter< real_t > streamWriter( mgr.vtkDataFormat_ );
      VTKMeshWriter::writePointsForMicroVertices( mgr.write2D_, streamWriter, storage, level );
      streamWriter.toStream( output );
   }

   output << "\n</DataArray>\n";
   output << "</Points>\n";

   if ( mgr.write2D_ )
   {
      VTKMeshWriter::writeCells2D( mgr.vtkDataFormat_, output, storage, levelinfo::num_microvertices_per_edge( level ) );
   }
   else
   {
      VTKMeshWriter::writeCells3D( mgr.vtkDataFormat_, output, storage, levelinfo::num_microvertices_per_edge( level ) );
   }

   output << "<CellData>\n";

   // write all scalar P0Functions of supported value type
   for ( const auto& function : mgr.feFunctionRegistry_.getP0Functions().getFunctions< double >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP0Functions().getFunctions< float >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP0Functions().getFunctions< int32_t >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getP0Functions().getFunctions< int64_t >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }

   output << "</CellData>\n";

   vtk::writePieceFooter( output );
}

template < typename value_t >
void VTKP0Writer::writeScalarFunction( std::ostream&                              output,
                                       const P0Function< value_t >&               function,
                                       const std::shared_ptr< PrimitiveStorage >& storage,
                                       const uint_t&                              level,
                                       bool                                       write2D,
                                       vtk::DataFormat                            vtkDataFormat )
{
   WALBERLA_ASSERT_EQUAL( storage, function.getStorage() );

   VTKStreamWriter< value_t > streamWriter( vtkDataFormat );
   vtk::openDataElement( output, typeToString< value_t >(), function.getFunctionName(), 1, vtkDataFormat );

   if ( write2D )
   {
      for ( auto& it : storage->getFaces() )
      {
         const auto  faceID = it.first;
         const auto& face   = *it.second;

         WALBERLA_CHECK_EQUAL( function.getDGFunction()->polynomialDegree( faceID ), 0 );
         WALBERLA_CHECK_EQUAL( function.getDGFunction()->basis()->numDoFsPerElement( 2, 0 ), 1 );

         const auto memLayout = function.getDGFunction()->volumeDoFFunction()->memoryLayout();
         auto       dofs      = function.getDGFunction()->volumeDoFFunction()->dofMemory( faceID, level );

         uint_t rowsize       = static_cast< uint_t >( levelinfo::num_microvertices_per_edge( level ) ) - 1u;
         uint_t inner_rowsize = rowsize;

         // NOTE: The order in which we export the cell values for the micro-faces needs to be
         //       identical to the order in which we enumerate them in VTKMeshWriter::writeConnectivityP1Triangles
         for ( uint_t i = 0; i < rowsize; ++i )
         {
            for ( uint_t j = 0; j < inner_rowsize - 1; ++j )
            {
               streamWriter << dofs[volumedofspace::indexing::index( j, i, facedof::FaceType::GRAY, 0, 1, level, memLayout )];
               streamWriter << dofs[volumedofspace::indexing::index( j, i, facedof::FaceType::BLUE, 0, 1, level, memLayout )];
            }
            streamWriter << dofs[volumedofspace::indexing::index( inner_rowsize - 1, i, facedof::FaceType::GRAY, 0, 1, level, memLayout )];
            --inner_rowsize;
         }
      }
   }
   else
   {
      for ( auto& it : storage->getCells() )
      {
         const auto  cellID = it.first;
         const auto& cell   = *it.second;

         WALBERLA_CHECK_EQUAL( function.getDGFunction()->polynomialDegree( cellID ), 0 );
         WALBERLA_CHECK_EQUAL( function.getDGFunction()->basis()->numDoFsPerElement( 3, 0 ), 1 );

         const auto memLayout = function.getDGFunction()->volumeDoFFunction()->memoryLayout();
         auto       dofs      = function.getDGFunction()->volumeDoFFunction()->dofMemory( cellID, level );

         for ( auto cellType : celldof::allCellTypes )
         {
            for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
            {
               streamWriter
                   << dofs[volumedofspace::indexing::index( idxIt.x(), idxIt.y(), idxIt.z(), cellType, 0, 1, level, memLayout )];
            }
         }
      }
   }

   streamWriter.toStream( output );

   output << "\n</DataArray>\n";
}

} // namespace hyteg
