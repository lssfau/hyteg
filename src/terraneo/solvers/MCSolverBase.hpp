/*
 * Copyright (c) 2025 Ponsuganth Ilangovan P
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

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {
template < typename BaseOperator_T >
class MCSolverBase
{
 public:
   MCSolverBase( const std::shared_ptr< PrimitiveStorage >& storage, const uint_t minLevel, const uint_t maxLevel )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   {}

   MCSolverBase( const std::shared_ptr< PrimitiveStorage >&        storage,
                 const uint_t                                      minLevel,
                 const uint_t                                      maxLevel,
                 const std::shared_ptr< Solver< BaseOperator_T > > solverPtr )
   : storage_( storage )
   , minLevel_( minLevel )
   , maxLevel_( maxLevel )
   , solverPtr_( solverPtr )
   {}

   std::shared_ptr< Solver< BaseOperator_T > > getSolver() { return solverPtr_; }

 protected:
   const std::shared_ptr< PrimitiveStorage >& storage_;

   const uint_t minLevel_;
   const uint_t maxLevel_;

   std::shared_ptr< Solver< BaseOperator_T > > solverPtr_;
};
} // namespace hyteg