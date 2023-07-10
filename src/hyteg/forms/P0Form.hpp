/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Marcus Mohr, Nils Kohl, Benjamin Mann.
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

class P0Form : public Form
{
 public:
   virtual ~P0Form() {}

   // 2D P0
   virtual void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const { WALBERLA_ABORT( "Not implemented." ); }
   virtual void integrateAll( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 1 >& elMat ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   }
   virtual void integrateRow( const uint_t& row, const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 1 >& elMat ) const;

   // 3D P0
   virtual void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const { WALBERLA_ABORT( "Not implemented." ); }
   virtual void integrateAll( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 1 >& elMat ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   }
   virtual void integrateRow( const uint_t& row, const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 1 >& elMat ) const;

 private:
   virtual void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 1 >& elMat ) const;
   virtual void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 1 >& elMat ) const;
};

} // namespace hyteg
