/*
 * Copyright (c) 2021 Benjamin Mann
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

// wrapper to enable combining P1 operators with higher order forms
template < class WrappedForm >
class P1Form : public WrappedForm
{
 public:
   void integrateRow( const uint_t& row, const std::array< Point3D, 3 >& coords, Matrixr< 1, 3 >& elMat ) const;
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

   void integrateRow( const uint_t& row, const std::array< Point3D, 4 >& coords, Matrixr< 1, 4 >& elMat ) const;
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
   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrixr< 1, 3 >& elMat ) const
   {
      Point3D row;
      this->integrate( coords, row );
      for ( int i = 0; i < 3; ++i )
      {
         elMat( 0, i ) = row[i];
      }
   }

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrixr< 1, 4 >& elMat ) const
   {
      Point4D row;
      this->integrate( coords, row );
      for ( int i = 0; i < 4; ++i )
      {
         elMat( 0, i ) = row[i];
      }
   }
};

} // namespace hyteg
