/*
* Copyright (c) 2024 Nils Kohl.
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

#include "hyteg/dataexport/VTKOutput/VTKPointCloudOutput.hpp"

#include "core/Format.hpp"
#include "core/mpi/MPITextFile.h"

#include "hyteg/dataexport/VTKOutput/VTKStreamWriter.hpp"

#include "vtk/UtilityFunctions.h"

namespace hyteg {

using walberla::vtk::typeToString;

void VTKPointCloudOutput::write( uint_t timestep ) const
{
   const std::string completeFilePath =
       walberla::format( "%s/%s%s.vtu", dir_.c_str(), filename_.c_str(), walberla::format( "_ts%u", timestep ).c_str() );

   std::ostringstream output;

   vtk::writeXMLHeader( output );

   // We are writing VTK cells - but those cells have type 'point' (as opposed to ,e.g., type 'triangle').
   // Thus, the number of cells and number of points are the same number.
   vtk::writePieceHeader( output, points_.size(), points_.size() );

   //////////////
   /// Points ///
   //////////////

   output << "<Points>\n";
   vtk::openDataElement( output, typeToString< real_t >(), "", 3, vtkDataFormat_ );

   {
      VTKStreamWriter< real_t > streamWriter( vtkDataFormat_ );

      for ( const auto& p : points_ )
      {
         streamWriter << p[0] << p[1] << p[2];
      }

      streamWriter.toStream( output );
   }

   output << "\n</DataArray>\n";
   output << "</Points>\n";

   /////////////
   /// Cells ///
   /////////////

   using CellType = uint32_t;

   output << "<Cells>\n";
   vtk::openDataElement( output, typeToString< CellType >(), "connectivity", 0, vtkDataFormat_ );

   {
      // connectivity
      VTKStreamWriter< CellType > streamWriterCells( vtkDataFormat_ );
      for ( CellType idx = 0; idx < points_.size(); idx++ )
      {
         streamWriterCells << idx;
      }
      streamWriterCells.toStream( output );
   }

   output << "\n</DataArray>\n";

   using OffsetType = uint32_t;

   vtk::openDataElement( output, typeToString< OffsetType >(), "offsets", 0, vtkDataFormat_ );

   {
      // offsets
      VTKStreamWriter< OffsetType > streamWriterOffsets( vtkDataFormat_ );
      for ( OffsetType idx = 0; idx < points_.size(); idx++ )
      {
         streamWriterOffsets << idx + OffsetType( 1 );
      }
      streamWriterOffsets.toStream( output );
   }

   output << "\n</DataArray>\n";

   using CellTypeType = uint16_t;

   vtk::openDataElement( output, typeToString< CellTypeType >(), "types", 0, vtkDataFormat_ );

   {
      // cell types
      VTKStreamWriter< CellTypeType > streamWriterTypes( vtkDataFormat_ );
      for ( size_t i = 0; i < points_.size(); ++i )
      {
         streamWriterTypes << 1;
      }
      streamWriterTypes.toStream( output );
   }

   output << "\n</DataArray>\n";
   output << "</Cells>\n";

   //////////////
   /// Values ///
   //////////////

   output << "<PointData>\n";

   for ( const auto& [key, data] : values_ )
   {
      VTKStreamWriter< real_t > streamWriter( vtkDataFormat_ );
      vtk::openDataElement( output, typeToString< real_t >(), key, 1, vtkDataFormat_ );

      for ( const auto& d : data )
      {
         streamWriter << d;
      }

      streamWriter.toStream( output );

      output << "\n</DataArray>\n";
   }

   output << "</PointData>\n";

   vtk::writePieceFooter( output );

   walberla::mpi::writeMPITextFile( completeFilePath, output.str() );

   WALBERLA_ROOT_SECTION()
   {
      std::ofstream pvtu_file;
      pvtu_file.open( completeFilePath.c_str(), std::ofstream::out | std::ofstream::app );
      WALBERLA_CHECK( !!pvtu_file, "[VTKWriter] Error opening file: " << completeFilePath );
      vtk::writeXMLFooter( pvtu_file );
      pvtu_file.close();
   }
}

} // namespace hyteg
