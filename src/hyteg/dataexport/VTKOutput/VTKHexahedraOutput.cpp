/*
* Copyright (c) 2025 Nils Kohl.
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

#include "VTKHexahedraOutput.hpp"

#include "core/mpi/MPITextFile.h"

#include "hyteg/Format.hpp"
#include "hyteg/types/PointND.hpp"

#include "VTKStreamWriter.hpp"
#include "vtk/UtilityFunctions.h"

namespace hyteg {

using walberla::uint_c;
using walberla::uint_t;

void VTKHexahedraOutput::addRectangle( const std::array< Point3D, 4 >& rectPoints )
{
   pointsRect_.push_back( rectPoints );
}

void VTKHexahedraOutput::addHexahedron( const std::array< Point3D, 8 >& hexahedronPoints )
{
   pointsHexa_.push_back( hexahedronPoints );
}

void VTKHexahedraOutput::write( uint_t timestep ) const
{
   const std::string completeFilePath =
       walberla::format( "%s/%s%s.vtu", dir_.c_str(), filename_.c_str(), walberla::format( "_ts%u", timestep ).c_str() );

   std::ostringstream output;

   vtk::writeXMLHeader( output );

   const auto numPoints = 4 * pointsRect_.size() + 8 * pointsHexa_.size();
   const auto numCells  = pointsRect_.size() + pointsHexa_.size();
   vtk::writePieceHeader( output, numPoints, numCells );

   //////////////
   /// Points ///
   //////////////

   output << "<Points>\n";
   vtk::openDataElement( output, walberla::vtk::typeToString< real_t >(), "", 3, vtkDataFormat_ );

   {
      VTKStreamWriter< real_t > streamWriter( vtkDataFormat_ );

      for ( const auto& rect : pointsRect_ )
      {
         for ( const auto& p : rect )
         {
            streamWriter << p[0] << p[1] << p[2];
         }
      }

      for ( const auto& hexa : pointsHexa_ )
      {
         for ( const auto& p : hexa )
         {
            streamWriter << p[0] << p[1] << p[2];
         }
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
   vtk::openDataElement( output, walberla::vtk::typeToString< CellType >(), "connectivity", 0, vtkDataFormat_ );

   {
      // connectivity
      VTKStreamWriter< CellType > streamWriterCells( vtkDataFormat_ );

      for ( CellType idx = 0; idx < numPoints; idx++ )
      {
         streamWriterCells << idx;
      }
      streamWriterCells.toStream( output );
   }

   output << "\n</DataArray>\n";

   using OffsetType = uint32_t;

   vtk::openDataElement( output, walberla::vtk::typeToString< OffsetType >(), "offsets", 0, vtkDataFormat_ );

   {
      // offsets
      VTKStreamWriter< OffsetType > streamWriterOffsets( vtkDataFormat_ );
      int                           offset = 0;
      for ( OffsetType idx = 0; idx < pointsRect_.size(); idx++ )
      {
         offset += 4;
         streamWriterOffsets << offset;
      }
      for ( OffsetType idx = 0; idx < pointsHexa_.size(); idx++ )
      {
         offset += 8;
         streamWriterOffsets << offset;
      }
      streamWriterOffsets.toStream( output );
   }

   output << "\n</DataArray>\n";

   using CellTypeType = uint16_t;

   vtk::openDataElement( output, walberla::vtk::typeToString< CellTypeType >(), "types", 0, vtkDataFormat_ );

   {
      // cell types
      VTKStreamWriter< CellTypeType > streamWriterTypes( vtkDataFormat_ );
      for ( OffsetType idx = 0; idx < pointsRect_.size(); idx++ )
      {
         streamWriterTypes << VTK_QUAD;
      }
      for ( OffsetType idx = 0; idx < pointsHexa_.size(); idx++ )
      {
         streamWriterTypes << VTK_HEXAHEDRON;
      }
      streamWriterTypes.toStream( output );
   }

   output << "\n</DataArray>\n";
   output << "</Cells>\n";

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
