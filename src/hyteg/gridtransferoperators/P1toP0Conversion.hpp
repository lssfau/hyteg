
/*
 * Copyright (c) 2025 Ponsuganth Ilangovan.
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

#include "hyteg/p0functionspace/P0Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/types/Averaging.hpp"

namespace hyteg {

/// Convert a P1 to a P0 function
///
/// The function takes a P2Function on level (P1Level-1) and converts it to a P1Function on level P1Level by
/// assigning the DoFs of the P2Function to those of the P1Function on the finer level.
void P1toP0Conversion( const P1Function< real_t >& src,
                       const P0Function< real_t >& dst,
                       const uint_t                P1Level,
                       const hyteg::AveragingType  averagingType );

} // namespace hyteg
