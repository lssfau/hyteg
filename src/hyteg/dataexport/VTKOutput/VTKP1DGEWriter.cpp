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
#include "hyteg/dataexport/VTKOutput/VTKP1DGEWriter.hpp"

#include "core/DataTypes.h"

#include "hyteg/dataexport/VTKOutput/VTKHelpers.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"

// from walberla
#include "vtk/UtilityFunctions.h"

namespace hyteg {

using walberla::vtk::typeToString;

// TODO:
// write all points (possibly also edges with option)
// write all cells (individually - double vertices - connect vertices from step 1)
// write all data via evaluation at the points (second DGFunction::evaluate( faceID, microIdx, microType ) would really speed up
// the process since we already know which element we want to evaluate - saves element search)

void VTKP1DGEWriter::write( const VTKOutput& mgr, std::ostream& output, const uint_t& level )
{
   if ( mgr.getNumRegisteredFunctions( vtk::DoFType::P1DGE ) == 0 )
   {
      return;
   }

   auto storage = mgr.storage_;

   const uint_t numberOfPoints2D = storage->getNumberOfLocalFaces() * levelinfo::num_microfaces_per_face( level ) * 3;
   const uint_t numberOfCells2D  = storage->getNumberOfLocalFaces() * levelinfo::num_microfaces_per_face( level );

   const uint_t numberOfPoints3D = storage->getNumberOfLocalCells() * levelinfo::num_microcells_per_cell( level ) * 4;
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

   VTKMeshWriter::writePointsForMicroVertices( mgr, output, storage, level, true );

   output << "\n</DataArray>\n";
   output << "</Points>\n";

   if ( mgr.write2D_ )
   {
      VTKMeshWriter::writeCells2D( mgr, output, storage, levelinfo::num_microvertices_per_edge( level ), true );
   }
   else
   {
      VTKMeshWriter::writeCells3D( mgr, output, storage, levelinfo::num_microvertices_per_edge( level ), true );
   }

   output << "<PointData>\n";

   // write all P1DGEVectorFunctions of supported value type
   for ( const auto& function : mgr.feFunctionRegistry_.getEGFunctions().getFunctions< double >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getEGFunctions().getFunctions< int32_t >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.feFunctionRegistry_.getEGFunctions().getFunctions< int64_t >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }

   output << "</PointData>\n";

   vtk::writePieceFooter( output );
}

template < typename value_t >
void VTKP1DGEWriter::writeVectorFunction( std::ostream&                              output,
                                          const EGFunction< value_t >&               function,
                                          const std::shared_ptr< PrimitiveStorage >& storage,
                                          const uint_t&                              level,
                                          bool                                       write2D,
                                          vtk::DataFormat                            vtkDataFormat )
{
   WALBERLA_ASSERT_EQUAL( storage, function.getStorage() );

   VTKOutput::VTKStreamWriter< value_t > streamWriter( vtkDataFormat );
   vtk::openDataElement( output, typeToString< value_t >(), function.getFunctionName(), function.getDimension(), vtkDataFormat );

   if ( write2D )
   {
      for ( const auto& it : storage->getFaces() )
      {
         const PrimitiveID faceID = it.first;
         const Face&       face   = *it.second;

         for ( auto faceType : facedof::allFaceTypes )
         {
            for ( const auto& idxIt : facedof::macroface::Iterator( level, faceType ) )
            {
               const std::array< indexing::Index, 3 > vertexIndices =
                   facedof::macroface::getMicroVerticesFromMicroFace( idxIt, faceType );
               for ( uint_t i = 0; i < 3; i++ )
               {
                  const auto vtkPoint = vertexdof::macroface::coordinateFromIndex( level, face, vertexIndices[i] );
                  value_t    value;
                  for ( uint_t j = 0; j < function.getDimension(); j += 1 )
                  {
                     function.evaluateOnMicroElement( vtkPoint, level, faceID, idxIt, faceType, j, value );
                     streamWriter << value;
                  }
               }
            }
         }
      }
   }
   else
   {
      for ( const auto& it : storage->getCells() )
      {
         const PrimitiveID cellID = it.first;
         const Cell&       cell   = *it.second;

         for ( auto cellType : celldof::allCellTypes )
         {
            for ( auto idxIt : celldof::macrocell::Iterator( level, cellType ) )
            {
               const std::array< indexing::Index, 4 > vertexIndices =
                   celldof::macrocell::getMicroVerticesFromMicroCell( idxIt, cellType );
               for ( uint_t i = 0; i < 4; i++ )
               {
                  const auto vtkPoint = vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[i] );
                  value_t    value;
                  for ( uint_t j = 0; j < 3; j += 1 )
                  {
                     function.evaluateOnMicroElement( vtkPoint, level, cellID, idxIt, cellType, j, value );
                     streamWriter << value;
                  }
               }
            }
         }
      }
   }

   streamWriter.toStream( output );

   output << "\n</DataArray>\n";
}

} // namespace hyteg
