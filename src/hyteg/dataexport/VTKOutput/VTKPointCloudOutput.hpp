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

#pragma once

#include "core/debug/CheckFunctions.h"

#include "hyteg/dataexport/VTKOutput/VTKHelpers.hpp"
#include "hyteg/types/PointND.hpp"

using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

/// \brief Simple class to write point clouds to VTK in parallel.
///
/// After setting a (process-local) point array, multiple value arrays can be added that are all written to VTK.
class VTKPointCloudOutput
{
 public:
   VTKPointCloudOutput( std::string dir, std::string filename )
   : dir_( dir )
   , filename_( filename )
   , vtkDataFormat_( vtk::DataFormat::BINARY ){};

   void setPoints( const std::vector< Point3D >& points ) { points_ = points; }

   void addPoints( const std::vector< Point3D >& points ) { points_.insert( points_.end(), points.begin(), points.end() ); }

   void setValues( std::string key, const std::vector< real_t >& values )
   {
      WALBERLA_CHECK_EQUAL( points_.size(), values.size(), "Values and point size arrays mismatch." )
      values_[key] = values;
   }

   void addValues( std::string key, const std::vector< real_t >& values )
   {
      values_[key].insert( values_[key].end(), values.begin(), values.end() );
   }

   void write( uint_t timestep = uint_c( 0 ) ) const;

 private:
   std::string dir_;
   std::string filename_;

   vtk::DataFormat vtkDataFormat_;

   std::vector< Point3D >                         points_;
   std::map< std::string, std::vector< real_t > > values_;
};

} // namespace hyteg