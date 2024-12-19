/*
 * Copyright (c) 2023 Daniel Bauer.
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

#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {
namespace n1e1 {

class N1E1Form_mass_LU_fp64
{
 public:
   /** Edges are in FEniCS ordering:
    * 
    *  3
    *  |\`\.
    *  | 0 `\.
    *  |  \   1
    *  3  2 _  `\.
    *  |  /  `-2 `\.
    *  | 4      `\_`\
    *  0------5------1
   **/
   void integrateAll( const std::array< PointND< double, 3 >, 4 >& coords,
                      const std::array< walberla::int16_t, 6 >&    edgeDirections,
                      Matrix< double, 6, 6 >&                      elMat ) const;
};

} // namespace n1e1
} // namespace hyteg
