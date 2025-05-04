/*
 * Copyright (c) 2025 Andreas Burkhart.
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
#include <map>

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/sparseassembly/VectorProxy.hpp"

namespace hyteg {

class MapVector : public VectorProxy
{
 public:
   MapVector() { valueMap = std::make_shared< std::map< uint_t, real_t > >(); }
   virtual ~MapVector() {}

   void setValue( uint_t idx, real_t value ) override { ( *valueMap )[idx] = value; }

   real_t getValue( uint_t idx ) const override { return ( *valueMap )[idx]; }

   uint_t getSize() const { return valueMap->size(); }

   std::shared_ptr< std::map< uint_t, real_t > > getValueMap() { return valueMap; }

 private:
   std::shared_ptr< std::map< uint_t, real_t > > valueMap;
};

} // namespace hyteg