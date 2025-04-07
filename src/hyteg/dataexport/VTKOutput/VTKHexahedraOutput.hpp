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

#pragma once

#include "core/debug/CheckFunctions.h"

#include "hyteg/dataexport/VTKOutput/VTKHelpers.hpp"
#include "hyteg/types/PointND.hpp"

using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

/// \brief Simple class to write rects and/or hexahedra to VTK.
class VTKHexahedraOutput
{
 public:
   VTKHexahedraOutput( std::string dir, std::string filename )
   : dir_( dir )
   , filename_( filename )
   , vtkDataFormat_( vtk::DataFormat::ASCII /* Somehow the binary output does not work for the hex output. */ ) {};

   /// Numbering is counter-clockwise.
   void addRectangle( const std::array< Point3D, 4 >& rectPoints );

   /// Numbering is counter-clockwise.
   /// For hexahedra: first "bottom" rect counter-clockwise, then "top" rect counter-clockwise.
   /// The starting point in the "top" layer must be the point "above" the starting point in the "bottom" layer.
   void addHexahedron( const std::array< Point3D, 8 >& hexahedronPoints );

   /// Writes for all processes. Duplicated rects and hexes will be written if duplicate rects
   /// and hexes have been added on many processes.
   ///
   /// Must be called collectively.
   void write( uint_t timestep = uint_c( 0 ) ) const;

 private:
   static constexpr int VTK_QUAD       = 9;  // For 2D rectangles (quads)
   static constexpr int VTK_HEXAHEDRON = 12; // For 3D hexahedra

   std::string dir_;
   std::string filename_;

   vtk::DataFormat vtkDataFormat_;

   std::vector< std::array< Point3D, 4 > > pointsRect_;
   std::vector< std::array< Point3D, 8 > > pointsHexa_;
};

} // namespace hyteg
