/*
 * Copyright (c) 2017-2021 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#pragma once

namespace hyteg {

namespace vtk {

enum class DataFormat
{
   ASCII,
   BINARY
};

enum class DoFType
{
   VERTEX,
   EDGE_X,
   EDGE_Y,
   EDGE_Z,
   EDGE_XY,
   EDGE_XZ,
   EDGE_YZ,
   EDGE_XYZ,
   DG,
   P2,
   N1E1
   P2,
   P1DGE
};

inline void writeXMLHeader( std::ostream& output )
{
   WALBERLA_ROOT_SECTION()
   {
      output << "<?xml version=\"1.0\"?>\n";
      output << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
      output << "<UnstructuredGrid>\n";
   }
}

inline void writeXMLFooter( std::ostream& output )
{
   WALBERLA_ROOT_SECTION()
   {
      output << "</UnstructuredGrid>\n";
      output << "</VTKFile>\n";
   }
}

inline void writePieceHeader( std::ostream& output, const uint_t& numberOfPoints, const uint_t& numberOfCells )
{
   output << "<Piece "
          << "NumberOfPoints=\"" << numberOfPoints << "\" "
          << "NumberOfCells=\"" << numberOfCells << "\""
          << ">\n";
}

inline void writePieceFooter( std::ostream& output )
{
   output << "</Piece>\n";
}

inline void
    openDataElement( std::ostream& output, const std::string& type, const std::string& name, uint_t nComponents, DataFormat fmt )
{
   // open element and write type
   output << "<DataArray type=\"" << type << "\"";

   // write name, if given
   if ( name.length() > 0 )
   {
      output << " Name=\"" << name << "\"";
   }

   // write number of components, if required
   if ( nComponents > 0 )
   {
      output << " NumberOfComponents=\"" << nComponents << "\"";
   }

   // specify format
   if ( fmt == DataFormat::ASCII )
   {
      output << " format=\"ascii\">\n";
   }
   else if ( fmt == DataFormat::BINARY )
   {
      output << " format=\"binary\">\n";
   }
   else
   {
      WALBERLA_ABORT( "VTK format not supported." );
   }
}

} // namespace vtk

} // namespace hyteg
