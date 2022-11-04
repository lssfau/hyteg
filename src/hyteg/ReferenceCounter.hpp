/*
* Copyright (c) 2022 Nils Kohl.
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
namespace internal {

using walberla::uint_t;

class ReferenceCounter
{
 public:
   ReferenceCounter()
   : refs_( 0 )
   {}

   void increaseRefs() { refs_++; }

   void decreaseRefs() { refs_--; }

   uint_t refs() const { return refs_; }

 private:
   uint_t refs_;
};

} // namespace internal
} // namespace hyteg