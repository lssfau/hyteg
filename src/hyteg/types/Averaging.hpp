/*
* Copyright (c) 2025 Ponsuganth Ilangovan
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

namespace hyteg {
///
/// \brief This enum flags can be used in P1toP0Conversion and P0toP0AveragedInjection 
///        to specify the averaging method that should to be used
///
enum class AveragingType
{
   ARITHMETIC, ///< Vertices and centroid (only vertices in 2D (P1 to P0 transfer), 
               ///< the first three flags can be used for transfer to P0 lower levels)
   HARMONIC,
   GEOMETRIC,
   ARITHMETIC_QP, ///< Suffix _QP referes to averaging on Quadrature points 
                  ///< (this and following is only applicable for P1 to P0 transfer)
   HARMONIC_QP,
   GEOMETRIC_QP
};
} // namespace hyteg