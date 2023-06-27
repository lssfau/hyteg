/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

namespace hyteg {

/// Base class for driver classes that export Finite Element functions
class FEFunctionWriter
{
 public:

   /// Add an FE Function to became part of the next dataexport phase
   template < template < typename > class func_t, typename value_t >
   inline void add( const func_t< value_t >& function )
   {
     WALBERLA_ABORT( "FEFunctionWriter::add() not overwritten in derived class!" );
   }

   /// Writes output only if writeFrequency > 0 and timestep % writeFrequency == 0.
   /// Therefore always writes output if timestep is 0.
   virtual void write( const uint_t& level, const uint_t& timestep = 0 ) const = 0;

   /// Set parameter specified by string key to the value specified by string value
   virtual void setParameter( const std::string& key, const std::string& value ) = 0;

};

} // namespace hyteg
