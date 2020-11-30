/*
 * Copyright (c) 2020 Marcus Mohr.
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

#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"

namespace hyteg {
namespace vectorFunctionTools {

/// From a vector of VectorFunctions extract a vector of their component functions for given dimension
template< typename vType >
std::vector< std::reference_wrapper< const typename vType::VectorComponentType > >
    filter( uint_t dim, const std::vector< std::reference_wrapper< const vType > >& functions )
{
   std::vector< std::reference_wrapper< const typename vType::VectorComponentType > > functions_scalar;
   std::transform( functions.begin(), functions.end(), std::back_inserter( functions_scalar ), [dim]( auto& function ) {
      return std::cref( function.get().component( dim ) );
   } );
   return functions_scalar;
}

} // namespace vectorFunctionTools
} // namespace hyteg
