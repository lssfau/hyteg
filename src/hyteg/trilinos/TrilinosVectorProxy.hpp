/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include "hyteg/sparseassembly/VectorProxy.hpp"
#include "hyteg/trilinos/TrilinosVector.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

using Teuchos::RCP;
using Teuchos::rcp;

template< typename TpetraVectorType >
class TrilinosVectorProxy : public VectorProxy
{
 public:

   explicit TrilinosVectorProxy( const RCP< TpetraVectorType >& vec )
   : vec_( vec )
   {}

   void setValue( uint_t idx, real_t value ) override { vec_->replaceGlobalValue( idx, 0, value ); }

   real_t getValue( uint_t idx ) const override
   {
      const auto             map        = vec_->getMap();
      const auto             localIndex = map->getLocalElement( idx );
      return vec_->getData( 0 ).operator[]( localIndex );
   }

 private:
   RCP< TpetraVectorType > vec_;
};

} // namespace hyteg