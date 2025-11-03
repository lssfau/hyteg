/*
* Copyright (c) 2025 Marcus Mohr.
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

#include <concepts>

#include "core/Environment.h"

#include "hyteg/functions/CSFVectorFunction.hpp"
#include "hyteg/functions/Function.hpp"

using walberla::real_t;

namespace hyteg::concepts {

/// Concept for admissable value types of finite element functions
template < typename T >
concept value_type = std::floating_point< T > || std::integral< T >;

/// Concept for macro-primitives
template < typename T >
concept concrete_primitive =
    std::is_same_v< T, Cell > || std::is_same_v< T, Face > || std::is_same_v< T, Edge > || std::is_same_v< T, Vertex >;

/// \concept fe_function_scalar
/// \brief Concept matching any kind of scalar-valued finite element function in HyTeG (T::getDimension() returns 1)
template < typename T >
concept fe_function_scalar =
    std::is_base_of_v< Function< T >, T > && !std::is_same_v< n1e1::N1E1VectorFunction< typename T::valueType >, T >;

/// Concept matching any kind of vector-valued finite element function in HyTeG
template < typename T >
concept fe_function_vectorial =
    std::is_base_of_v< CSFVectorFunction< T >, T > || std::is_same_v< n1e1::N1E1VectorFunction< typename T::valueType >, T >;

/// Concept matching any kind of block/composite finite element function in HyTeG (e.g. P2-P1 Taylor-Hood)
template < typename T >
concept fe_function_composite = std::is_base_of_v< BlockFunction< typename T::ValueType >, T > ||
                                std::is_base_of_v< BlockFunction< typename T::valueType >, T >;

/// Concept matching any kind of Finite Element function in HyTeG
template < typename T >
concept fe_function = fe_function_scalar< T > || fe_function_vectorial< T > || fe_function_composite< T >;

} // namespace hyteg::concepts
