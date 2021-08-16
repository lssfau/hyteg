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
#include "hyteg/dataexport/VTKDGDoFWriter.hpp"

#include "core/DataTypes.h"

#include "hyteg/dataexport/VTKHelpers.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"

// from walberla
#include "vtk/UtilityFunctions.h"

namespace hyteg {

using walberla::vtk::typeToString;

void VTKDGDoFWriter::write( const VTKOutput& mgr, std::ostream& output, const uint_t& level )
{
   if ( mgr.dgFunctions_.size() == 0 )
   {
      return;
   }

   auto storage = mgr.storage_;

   const uint_t numberOfPoints = storage->getNumberOfLocalFaces() * levelinfo::num_microvertices_per_face( level );
   const uint_t numberOfCells  = storage->getNumberOfLocalFaces() * levelinfo::num_microfaces_per_face( level );

   vtk::writePieceHeader( output, numberOfPoints, numberOfCells );

   output << "<Points>\n";
   vtk::openDataElement( output, typeToString< real_t >(), "", 3, mgr.vtkDataFormat_ );

   VTKMeshWriter::writePointsForMicroVertices( mgr, output, storage, level );

   output << "\n</DataArray>\n";
   output << "</Points>\n";

   VTKMeshWriter::writeCells2D( mgr, output, storage, levelinfo::num_microvertices_per_edge( level ) );

   output << "<CellData>";

   for ( const auto& function : mgr.dgFunctions_.getFunctions< double >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.dgFunctions_.getFunctions< int32_t >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.dgFunctions_.getFunctions< int64_t >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }

   output << "\n</CellData>\n";

   vtk::writePieceFooter( output );
}

template < typename value_t >
void VTKDGDoFWriter::writeScalarFunction( std::ostream&                              output,
                                          const DGFunction< value_t >&               function,
                                          const std::shared_ptr< PrimitiveStorage >& storage,
                                          const uint_t&                              level,
                                          bool                                       write2D,
                                          vtk::DataFormat                            vtkDataFormat )
{
   vtk::openDataElement( output, typeToString< real_t >(), function.getFunctionName(), 1, vtkDataFormat );

   for ( const auto& it : storage->getFaces() )
   {
      const Face& face = *it.second;

      uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );
      uint_t inner_rowsize = rowsize;
      output << std::scientific;

      uint_t idx;

      for ( size_t j = 0; j < rowsize - 1; ++j )
      {
         for ( size_t i = 0; i < inner_rowsize - 2; ++i )
         {
            idx = facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_GRAY_C );
            output << face.getData( function.getFaceDataID() )->getPointer( level )[idx] << " ";
            idx = facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_BLUE_C );
            output << face.getData( function.getFaceDataID() )->getPointer( level )[idx] << " ";
         }
         idx = facedof::macroface::indexFaceFromGrayFace( level, inner_rowsize - 2, j, stencilDirection::CELL_GRAY_C );
         output << face.getData( function.getFaceDataID() )->getPointer( level )[idx] << " ";
         --inner_rowsize;
      }
   }
   output << "\n</DataArray>\n";
}

} // namespace hyteg
