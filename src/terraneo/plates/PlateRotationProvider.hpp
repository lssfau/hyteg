/*
 * Copyright (c) 2022 Berta Vilacis, Marcus Mohr.
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

#include "terraneo/dataimport/FileIO.hpp"

namespace terraneo {
namespace plates {

/// Class storing and providing access to plate rotation information
///
/// This class deals with importating infromation reconstructed plate
/// rotation information from an external datafile, stores that information
/// and provides access to it.
class PlateRotationProvider
{
 public:
   template < typename ImportStrategy >
   PlateRotationProvider( std::string nameOfRotationsFile, ImportStrategy readRotationsFile )
   {
      rotInfos_ = readRotationsFile( nameOfRotationsFile );
   }

   const std::vector< RotationInfo >& getRotations() const { return rotInfos_; };

 private:
   std::vector< RotationInfo > rotInfos_;
};

} // namespace plates
} // namespace terraneo
