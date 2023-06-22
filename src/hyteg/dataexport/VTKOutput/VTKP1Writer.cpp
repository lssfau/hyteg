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
#include "hyteg/dataexport/VTKOutput/VTKP1Writer.hpp"

#include "core/DataTypes.h"

#include "hyteg/dataexport/VTKOutput/VTKHelpers.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"

// from walberla
#include "vtk/UtilityFunctions.h"

namespace hyteg {

using walberla::vtk::typeToString;

void VTKP1Writer::write( const VTKOutput& mgr, std::ostream& output, const uint_t& level )
{
   if ( mgr.getNumRegisteredFunctions( vtk::DoFType::VERTEX ) == 0 )
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

   VTKMeshWriter::writePointsForMicroVertices( mgr, output, storage, level );

   output << "\n</DataArray>\n";
   output << "</Points>\n";

   if ( mgr.write2D_ )
   {
      VTKMeshWriter::writeCells2D( mgr, output, storage, levelinfo::num_microvertices_per_edge( level ) );
   }
   else
   {
      VTKMeshWriter::writeCells3D( mgr, output, storage, levelinfo::num_microvertices_per_edge( level ) );
   }

   output << "<PointData>\n";

   // write all scalar P1Functions of supported value type
   for ( const auto& function : mgr.p1Functions_.getFunctions< double >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.p1Functions_.getFunctions< int32_t >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.p1Functions_.getFunctions< int64_t >() )
   {
      writeScalarFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }

   for ( const auto& function : mgr.p1VecFunctions_.getFunctions< double >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.p1VecFunctions_.getFunctions< int32_t >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }
   for ( const auto& function : mgr.p1VecFunctions_.getFunctions< int64_t >() )
   {
      writeVectorFunction( output, function, storage, level, mgr.write2D_, mgr.vtkDataFormat_ );
   }

   output << "</PointData>\n";

   vtk::writePieceFooter( output );
}

template < typename value_t >
void VTKP1Writer::writeScalarFunction( std::ostream&                                  output,
                                       const vertexdof::VertexDoFFunction< value_t >& function,
                                       const std::shared_ptr< PrimitiveStorage >&     storage,
                                       const uint_t&                                  level,
                                       bool                                           write2D,
                                       vtk::DataFormat                                vtkDataFormat )
{
   WALBERLA_ASSERT_EQUAL( storage, function.getStorage() );

   VTKOutput::VTKStreamWriter< value_t > streamWriter( vtkDataFormat );
   vtk::openDataElement( output, typeToString< value_t >(), function.getFunctionName(), 1, vtkDataFormat );

   if ( write2D )
   {
      for ( const auto& it : storage->getFaces() )
      {
         const Face& face = *it.second;

         size_t len = levelinfo::num_microvertices_per_face( level );

         for ( size_t i = 0; i < len; ++i )
         {
            streamWriter << face.getData( function.getFaceDataID() )->getPointer( level )[i];
         }
      }
   }
   else
   {
      for ( const auto& it : storage->getCells() )
      {
         const Cell& cell     = *it.second;
         const auto  cellData = cell.getData( function.getCellDataID() )->getPointer( level );

         for ( const auto& idxIt : vertexdof::macrocell::Iterator( level ) )
         {
            streamWriter << cellData[vertexdof::macrocell::index( level, idxIt.x(), idxIt.y(), idxIt.z() )];
         }
      }
   }

   streamWriter.toStream( output );

   output << "\n</DataArray>\n";
}

template < typename value_t >
void VTKP1Writer::writeVectorFunction( std::ostream&                              output,
                                       const P1VectorFunction< value_t >&         function,
                                       const std::shared_ptr< PrimitiveStorage >& storage,
                                       const uint_t&                              level,
                                       bool                                       write2D,
                                       vtk::DataFormat                            vtkDataFormat )
{
   WALBERLA_ASSERT_EQUAL( storage, function.getStorage() );

   VTKOutput::VTKStreamWriter< value_t > streamWriter( vtkDataFormat );

   uint_t dim = function.getDimension();
   vtk::openDataElement( output, typeToString< value_t >(), function.getFunctionName(), dim, vtkDataFormat );

   if ( write2D )
   {
      for ( const auto& it : storage->getFaces() )
      {
         const Face& face = *it.second;

         size_t len = levelinfo::num_microvertices_per_face( level );

         for ( size_t i = 0; i < len; ++i )
         {
            for ( uint_t idx = 0; idx < dim; ++idx )
            {
               streamWriter << face.getData( function[idx].getFaceDataID() )->getPointer( level )[i];
            }
            // streamWriter << real_t(0); Paraview needs 3D vector fields to form glyphs
         }
      }
   }
   else
   {
      for ( const auto& it : storage->getCells() )
      {
         const Cell& cell      = *it.second;
         const auto  cellData0 = cell.getData( function[0].getCellDataID() )->getPointer( level );
         const auto  cellData1 = cell.getData( function[1].getCellDataID() )->getPointer( level );
         const auto  cellData2 = cell.getData( function[2].getCellDataID() )->getPointer( level );

         for ( const auto& idxIt : vertexdof::macrocell::Iterator( level ) )
         {
            streamWriter << cellData0[vertexdof::macrocell::index( level, idxIt.x(), idxIt.y(), idxIt.z() )];
            streamWriter << cellData1[vertexdof::macrocell::index( level, idxIt.x(), idxIt.y(), idxIt.z() )];
            streamWriter << cellData2[vertexdof::macrocell::index( level, idxIt.x(), idxIt.y(), idxIt.z() )];
         }
      }
   }

   streamWriter.toStream( output );

   output << "\n</DataArray>\n";
}

} // namespace hyteg
