/*
 * Copyright (c) 2023 Marcus Mohr.
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

using walberla::uint_t;

/// Base class for driver classes that store checkpoints
template < typename exporter_t >
class CheckpointExporter
{
 public:
   /// Register an FE Function to be included into checkpoints
   ///
   /// By calling this method the passed function object will be included into all future checkpoints.
   /// Data will be stored for all levels from minLevel up to maxLevel.
   template < template < typename > class func_t, typename value_t >
   inline void registerFunction( const func_t< value_t >& function, uint_t minLevel, uint_t maxLevel )
   {
      static_cast< exporter_t* >( this )->registerFunction( function, minLevel, maxLevel );
   }

#ifdef FE_FUNCTION_REGISTRY_HAS_REMOVE
   /// Deregister an FE Function to be no longer included into checkpoints
   ///
   /// By calling this method the passed function object will be excluded from future checkpoints.
   template < template < typename > class func_t, typename value_t >
   inline void deregisterFunction( const func_t< value_t >& function )
   {
      static_cast< exporter_t* >( this )->deregisterFunction( function );
   }
#endif

   /// Trigger storing of a single checkpoint
   ///
   /// \param filePath      Path to directory where the BP files are stored
   /// \param fileName      Name for output file
   inline void storeCheckpoint( std::string filePath, std::string fileName )
   {
      static_cast< exporter_t* >( this )->storeCheckpoint( filePath, fileName );
   }
};

} // namespace hyteg
