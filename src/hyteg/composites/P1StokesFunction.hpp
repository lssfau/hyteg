/*
 * Copyright (c) 2017-2022 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/functions/BlockFunction.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"

namespace hyteg {

template < typename ValueType >
class P1StokesFunction : public BlockFunction< ValueType >
{
 public:
   using valueType = ValueType;

   template < typename VType >
   using FunctionType = P1StokesFunction< VType >;

   using VelocityFunction_T = P1VectorFunction< ValueType >;
   using PressureFunction_T = P1Function< ValueType >;

   using Tag = typename FunctionTrait< P1StokesFunction< ValueType > >::Tag;

   P1StokesFunction( const std::string&                         _name,
                     const std::shared_ptr< PrimitiveStorage >& storage,
                     size_t                                     minLevel,
                     size_t                                     maxLevel )
   : P1StokesFunction( _name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
   {}

   P1StokesFunction( const std::string&                         _name,
                     const std::shared_ptr< PrimitiveStorage >& storage,
                     size_t                                     minLevel,
                     size_t                                     maxLevel,
                     BoundaryCondition                          velocityBC )
   : BlockFunction< ValueType >( _name )
   {
      this->subFunc_.push_back( std::make_shared< FunctionWrapper< P1VectorFunction< ValueType > > >(
          _name + "_uvw", storage, minLevel, maxLevel, velocityBC ) );
      this->subFunc_.push_back( std::make_shared< FunctionWrapper< P1Function< ValueType > > >(
          _name + "_p", storage, minLevel, maxLevel, BoundaryCondition::createAllInnerBC() ) );
   }

   [[nodiscard]] const P1VectorFunction< ValueType >& uvw() const {
      return this->subFunc_[0]->template unwrap< P1VectorFunction< ValueType > >();
   }

   [[nodiscard]] P1VectorFunction< ValueType >& uvw()
   {
      return this->subFunc_[0]->template unwrap< P1VectorFunction< ValueType > >();
   }

   [[nodiscard]] const P1Function< ValueType >& p() const {
      return this->subFunc_[1]->template unwrap< P1Function< ValueType > >();
   }

   [[nodiscard]] P1Function< ValueType >& p()
   {
      return this->subFunc_[1]->template unwrap< P1Function< ValueType > >();
   }

   /// \todo Get rid of this
   bool isDummy() const { return false; }
};

inline unsigned long long p1p1localFunctionMemorySize( const uint_t& level, const std::shared_ptr< PrimitiveStorage >& storage )
{
   if ( storage->hasGlobalCells() )
   {
      return 4 * vertexDoFLocalFunctionMemorySize( level, storage );
   }
   else
   {
      return 3 * vertexDoFLocalFunctionMemorySize( level, storage );
   }
}

inline unsigned long long p1p1globalFunctionMemorySize( const uint_t& level, const std::shared_ptr< PrimitiveStorage >& storage )
{
   const auto memLocal  = p1p1localFunctionMemorySize( level, storage );
   const auto memGlobal = walberla::mpi::allReduce( memLocal, walberla::mpi::SUM );
   return memGlobal;
}

} // namespace hyteg
