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
#pragma once

#include "core/Abort.h"

#include "hyteg/forms/Form.hpp"

namespace hyteg {

class P1Form : public Form
{
 public:
   virtual ~P1Form() {}

   // 2D P1
   virtual void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const { WALBERLA_ABORT( "Not implemented." ); }
   virtual void integrateAll( const std::array< Point3D, 3 >& coords, Matrix3r& elMat ) const { WALBERLA_ABORT( "Not implemented." ); }
   virtual void integrateRow( const uint_t & row, const std::array< Point3D, 3 >& coords, Matrixr< 1, 3 >& elMat ) const;

   // 3D P1
   virtual void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const { WALBERLA_ABORT( "Not implemented." ); }
   virtual void integrateAll( const std::array< Point3D, 4 >& coords, Matrix4r& elMat ) const { WALBERLA_ABORT( "Not implemented." ); }
   virtual void integrateRow( const uint_t & row, const std::array< Point3D, 4 >& coords, Matrixr< 1, 4 >& elMat ) const;

 private:

   virtual void integrateRow0( const std::array< Point3D, 3 >& coords, Matrixr< 1, 3 >& elMat ) const;
   virtual void integrateRow0( const std::array< Point3D, 4 >& coords, Matrixr< 1, 4 >& elMat ) const;

};

} // namespace hyteg
