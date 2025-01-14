/*
 * Copyright (c) 2025 Benjamin Mann.
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

#include "polynomial.hpp"

namespace hyteg {
namespace surrogate {
namespace polynomial {

/** Initialization of bases
 *  todo: add bases for higher degree spaces if needed
 */
const std::vector< Basis > Polynomial::basis_q = { Basis( 0 ),
                                                   Basis( 1 ),
                                                   Basis( 2 ),
                                                   Basis( 3 ),
                                                   Basis( 4 ),
                                                   Basis( 5 ),
                                                   Basis( 6 ),
                                                   Basis( 7 ),
                                                   Basis( 8 ),
                                                   Basis( 9 ),
                                                   Basis( 10 ),
                                                   Basis( 11 ),
                                                   Basis( 12 ) };

} // namespace polynomial
} // namespace surrogate
} // namespace hyteg
