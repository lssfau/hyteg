/*
 * Copyright (c) 2017-2020 Marcus Mohr.
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
#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

class DG1ToP2PlusBubbleFormHyTeG : public Form
{
 public:
   virtual ~DG1ToP2PlusBubbleFormHyTeG() {}

   virtual void integrateAll( const std::array< Point3D, 3 >& coords, Matrixr< 7, 3 >& elMat ) const = 0;

   virtual void integrateAll( const std::array< Point3D, 4 >& coords, Matrixr< 15, 4 >& elMat ) const
   {
      WALBERLA_ABORT( "not implemented" );
   };

 private:
   virtual void integrateRow0( const std::array< Point3D, 3 >& coords, Matrix< real_t, 1, 3 >& elMat ) const
   {
      WALBERLA_ABORT( "not implemented" );
   };

   virtual void integrateRow0( const std::array< Point3D, 4 >& coords, Matrix< real_t, 1, 4 >& elMat ) const
   {
      WALBERLA_ABORT( "not implemented" );
   };
};

} // namespace hyteg
