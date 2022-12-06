/*
 * Copyright (c) 2017-2019 Marcus Mohr.
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

#include "hyteg/misc/zeros.hpp"

#include <core/DataTypes.h>
#include <core/logging/Logging.h>

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

template < typename ValueType >
ValueType generateZero()
{
   WALBERLA_ABORT( "Don't know how to generate a 'zero' for given ValueType!" );
}

template <>
double generateZero< double >()
{
   return 0.0;
}

template <>
float generateZero< float >()
{
   return 0.0f;
}

template <>
uint_t generateZero< uint_t >()
{
   return 0;
}

template <>
int generateZero< int >()
{
   return 0;
}

template <>
long generateZero< long >()
{
   return 0;
}

template <>
long long generateZero< long long >()
{
   return 0;
}

} // namespace hyteg
