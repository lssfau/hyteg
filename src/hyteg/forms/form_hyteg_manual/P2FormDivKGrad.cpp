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
 * along with this progra m. If not, see <http://www.gnu.org/licenses/>.
 */

#include "hyteg/forms/form_hyteg_manual/P2FormDivKGrad.hpp"

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

real_t neutralCallback( const Point3D& x )
{
   WALBERLA_UNUSED( x );
   return real_c( 1.0 );
}

std::function< real_t( const Point3D& ) > P2Form_divKgrad::callback = &neutralCallback;

} // namespace hyteg
