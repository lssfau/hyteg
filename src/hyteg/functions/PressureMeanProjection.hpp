/*
* Copyright (c) 2024 Andreas Burkhart.
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

#include "core/DataTypes.h"

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/types/types.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

template < class SrcType >
inline void projectPressureMean( const SrcType& pressure, const uint_t& level )
{}

template <>
inline void projectPressureMean( const hyteg::P1Function< real_t >& pressure, const uint_t& level )
{
   hyteg::vertexdof::projectMean( pressure, level );
}

// template specialisations for additional function types go here

} // namespace hyteg