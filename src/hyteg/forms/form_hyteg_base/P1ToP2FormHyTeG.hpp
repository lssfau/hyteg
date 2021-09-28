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

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/types/matrix.hpp"
#include "hyteg/types/pointnd.hpp"
#include "hyteg/forms/Form.hpp"

namespace hyteg {

class P1ToP2FormHyTeG : public Form
{
 public:
   virtual ~P1ToP2FormHyTeG() {}

   virtual void integrateAll( const std::array< Point3D, 3 >& coords, Matrixr< 6, 3 >& elMat ) const = 0;

   virtual void integrateAll( const std::array< Point3D, 4 >& coords, Matrixr< 10, 4 >& elMat ) const = 0;
};

} // namespace hyteg
