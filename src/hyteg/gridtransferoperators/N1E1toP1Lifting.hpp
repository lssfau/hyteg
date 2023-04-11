/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include "hyteg/n1e1functionspace/N1E1VectorFunction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"

namespace hyteg {
namespace n1e1 {

/// Determine the scalar potential of an N1E1VectorFunction as a ::P1Function.
/// This is the adjoint of the gradient and equal to \f$-\mathrm{div}\f$.
void N1E1toP1Lifting( const N1E1VectorFunction< real_t >& src,
                      const P1Function< real_t >&         dst,
                      const uint_t                        lvl,
                      const DoFType&                      flag = All );

} // namespace n1e1
} // namespace hyteg
