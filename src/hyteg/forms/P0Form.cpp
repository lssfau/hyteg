/*
 * Copyright (c) 2017-2021 Dominik Thoennes, Marcus Mohr, Nils Kohl, Benjamin Mann.
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

#include "hyteg/forms/P0Form.hpp"

namespace hyteg {

void P0Form::integrateRow( const uint_t& row, const std::array< Point3D, 3 >& coords, Matrixr< 1, 1 >& elMat ) const
{
   switch ( row )
   {
   case 0:
      integrateRow0( coords, elMat );
      break;
   default:
      WALBERLA_ABORT( "P0Form::integrateRow() not implemented for row " << row );
   }
}

void P0Form::integrateRow( const uint_t& row, const std::array< Point3D, 4 >& coords, Matrixr< 1, 1 >& elMat ) const
{
   switch ( row )
   {
   case 0:
      integrateRow0( coords, elMat );
      break;
   default:
      WALBERLA_ABORT( "P0Form::integrateRow() not implemented for row " << row );
   }
}

void P0Form::integrateRow0( const std::array< Point3D, 3 >& coords, Matrixr< 1, 1 >& elMat ) const
{
   WALBERLA_UNUSED( coords );
   WALBERLA_UNUSED( elMat );
   WALBERLA_ABORT( "P0Form::integrateRow() not implemented for row 0" );
}

void P0Form::integrateRow0( const std::array< Point3D, 4 >& coords, Matrixr< 1, 1 >& elMat ) const
{
   WALBERLA_UNUSED( coords );
   WALBERLA_UNUSED( elMat );
   WALBERLA_ABORT( "P0Form::integrateRow() not implemented for row 0" );
}

} // namespace hyteg
