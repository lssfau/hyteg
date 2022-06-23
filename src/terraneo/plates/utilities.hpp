/*
 * Copyright (c) 2022 Marcus Mohr.
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

namespace terraneo {
namespace plates {

/// Convert a key string (stage name) to a plate stage (age value)
inline real_t keyStrToAge( const std::string& key )
{
   std::string ageStr = key.substr( 9, 6 );
   return PLATES_IO_STR_TO_FP( ageStr );
}

/// Convert a plate stage (age value) to a key string (stage name)
inline std::string ageToKeyStr( const real_t age )
{
   std::stringstream key;
   key.precision( 4 );
   key << std::fixed;
   key << "topology_" << age << "Ma_polygon";
   return key.str();
}

} // namespace plates
} // namespace terraneo
