/*
 * Copyright (c) 2021 Benjamin Mann
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

#include "core/Abort.h"
#include "core/DataTypes.h"

#include "hyteg/forms/P1Form.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

// wrapper to enable combining P1 operators with higher order forms
template < class WrappedForm >
class P1WrapperForm : public Form
{
 public:
   // use Ctor from WrappedForm
   using Form::Form;

   // conversion WrappedForm -> P1Wrapper
   P1WrapperForm( const WrappedForm& wf )
   : wrappedForm_( wf )
   {}
   // conversion WrappedForm -> P1Wrapper
   P1WrapperForm( WrappedForm&& wf )
   : wrappedForm_( wf )
   {}

   //   template < typename T = WrappedForm >
   //   typename std::enable_if< std::is_same< T, P1Form >::value, void >::type
   void integrateRow( const walberla::uint_t& row, const std::array< Point3D, 3 >& coords, Matrixr< 1, 3 >& elMat ) const
   {
      switch ( row )
      {
      case 0:
         integrateRow0( coords, elMat );
         break;
      default:
         WALBERLA_ABORT( "P1Form::integrateRow() not implemented for row " << row );
      }
   }

   //   template < typename T = WrappedForm >
   //   typename std::enable_if< std::is_same< T, P1Form >::value, void >::type
   void integrateRow( const walberla::uint_t& row, const std::array< Point3D, 4 >& coords, Matrixr< 1, 4 >& elMat ) const
   {
      switch ( row )
      {
      case 0:
         integrateRow0( coords, elMat );
         break;
      default:
         WALBERLA_ABORT( "P1Form::integrateRow() not implemented for row " << row );
      }
   }

 private:
   WrappedForm wrappedForm_;

   // extract vertex to vertex DoF part of WrappedForm
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrixr< 1, 3 >& elMat ) const
   {
      Point3D row;
      wrappedForm_.integrate( coords, row );
      for ( int i = 0; i < 3; ++i )
      {
         elMat( 0, i ) = row[i];
      }
   }

   // extract vertex to vertex DoF part of WrappedForm
   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrixr< 1, 4 >& elMat ) const
   {
      Point4D row;
      wrappedForm_.integrate( coords, row );
      for ( int i = 0; i < 4; ++i )
      {
         elMat( 0, i ) = row[i];
      }
   }
};

} // namespace hyteg
