/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include "hyteg/dataexport/VTKOutput/VTKN1E1Writer.hpp"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/n1e1functionspace/N1E1Indexing.hpp"
#include "hyteg/volumedofspace/CellDoFIndexing.hpp"

#include "vtk/UtilityFunctions.h"

namespace hyteg {

using walberla::vtk::typeToString;

void VTKN1E1Writer::write( const VTKOutput& mgr, std::ostream& output, const uint_t& level )
{
   if ( mgr.getNumRegisteredFunctions( vtk::DoFType::N1E1 ) == 0 )
   {
      return;
   }
   if ( mgr.write2D_ )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Skipping 3D N1E1 vector functions" )
      return;
   }

   auto storage = mgr.storage_;

   const uint_t numberOfPoints = storage->getNumberOfLocalCells() * levelinfo::num_microvertices_per_cell( level );
   const uint_t numberOfCells  = storage->getNumberOfLocalCells() * levelinfo::num_microcells_per_cell( level );

   vtk::writePieceHeader( output, numberOfPoints, numberOfCells );

   output << "<Points>\n";
   vtk::openDataElement( output, typeToString< real_t >(), "", 3, mgr.vtkDataFormat_ );

   VTKMeshWriter::writePointsForMicroVertices( mgr, output, storage, level );

   output << "\n</DataArray>\n";
   output << "</Points>\n";

   VTKMeshWriter::writeCells3D( mgr, output, storage, levelinfo::num_microvertices_per_edge( level ) );

   output << "<CellData>\n";

   // write all N1E1VectorFunctions of supported value type
   for ( const auto& function : mgr.n1e1Functions_.getFunctions< double >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.n1e1Functions_.getFunctions< int32_t >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.n1e1Functions_.getFunctions< int64_t >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.vtkDataFormat_ );
   }

   output << "</CellData>\n";

   vtk::writePieceFooter( output );
}

template < typename value_t >
void VTKN1E1Writer::writeVectorFunction( std::ostream&                              output,
                                         const n1e1::N1E1VectorFunction< value_t >& function,
                                         const std::shared_ptr< PrimitiveStorage >& storage,
                                         const uint_t&                              level,
                                         vtk::DataFormat                            vtkDataFormat )
{
   WALBERLA_ASSERT_EQUAL( storage, function.getStorage() );

   vtk::openDataElement( output, typeToString< value_t >(), function.getFunctionName(), function.getDimension(), vtkDataFormat );

   VTKOutput::VTKStreamWriter< value_t > streamWriter( vtkDataFormat );

   for ( const auto& it : storage->getCells() )
   {
      const PrimitiveID cellId = it.first;
      const Cell&       cell   = *it.second;

      for ( auto cellType : celldof::allCellTypes )
      {
         for ( const auto& idxIt : celldof::macrocell::Iterator( level, cellType ) )
         {
            const std::array< indexing::Index, 4 > vertexIndices =
                n1e1::macrocell::getMicroVerticesFromMicroCell( idxIt, cellType );
            const Point3D vtkPoint = 0.25 * ( vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[0] ) +
                                              vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[1] ) +
                                              vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[2] ) +
                                              vertexdof::macrocell::coordinateFromIndex( level, cell, vertexIndices[3] ) );

            typename n1e1::N1E1VectorFunction< value_t >::VectorType value;
            function.evaluateOnMicroElement( vtkPoint, level, cellId, idxIt, cellType, value );

            streamWriter << value[0];
            streamWriter << value[1];
            streamWriter << value[2];
         }
      }
   }

   streamWriter.toStream( output );

   output << "\n</DataArray>\n";
}

} // namespace hyteg
