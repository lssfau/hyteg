/*
* Copyright (c) 2017-2024 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/dataexport/VTKOutput/VTKHelpers.hpp"

#include "vtk/Base64Writer.h"

namespace hyteg {

/// Wrapper class that handles writing data in ASCII or binary format.
///
/// \tparam DTypeInVTK data type that the input data is converted to before writing it to the VTK file
template < typename DTypeInVTK >
class VTKStreamWriter
{
 public:
   explicit VTKStreamWriter( vtk::DataFormat vtkDataFormat )
   : vtkDataFormat_( vtkDataFormat )
   {
      if ( vtkDataFormat_ == vtk::DataFormat::ASCII )
      {
         outputAscii_ << std::scientific;
      }
   }

   template < typename T >
   VTKStreamWriter& operator<<( const T& data )
   {
      if ( vtkDataFormat_ == vtk::DataFormat::ASCII )
      {
         outputAscii_ << static_cast< DTypeInVTK >( data ) << "\n";
      }
      else if ( vtkDataFormat_ == vtk::DataFormat::BINARY )
      {
         outputBase64_ << static_cast< DTypeInVTK >( data );
      }

      return *this;
   }

   void toStream( std::ostream& os )
   {
      if ( vtkDataFormat_ == vtk::DataFormat::ASCII )
      {
         os << outputAscii_.str();
         // reset string stream
         // outputAscii_.str( std::string() );
         // outputAscii_.clear();
      }
      else if ( vtkDataFormat_ == vtk::DataFormat::BINARY )
      {
         outputBase64_.toStream( os );
         // Base64Writer::toStream() already reset the object
         // so nothing left to do for us here
      }
   }

 private:
   vtk::DataFormat             vtkDataFormat_;
   std::ostringstream          outputAscii_;
   walberla::vtk::Base64Writer outputBase64_;
};

} // namespace hyteg
