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

#include "hyteg/forms/form_hyteg_generated/P2FormDivKGradBlending.hpp"

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/types/matrix.hpp"
#include "hyteg/types/pointnd.hpp"

namespace hyteg {

real_t neutralCallback( const Point3D& x, const std::shared_ptr< GeometryMap >& map )
{
   WALBERLA_UNUSED( x );
   WALBERLA_UNUSED( map );
   return real_c( 1.0 );
}

std::function< real_t( const Point3D&, const std::shared_ptr< GeometryMap >& ) > P2Form_divKgradBlending::callback = &neutralCallback;

} // namespace hyteg
