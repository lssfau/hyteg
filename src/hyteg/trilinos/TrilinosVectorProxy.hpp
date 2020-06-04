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

class TrilinosVectorProxy : public VectorProxy
{
 public:
   TrilinosVectorProxy( const RCP< Tpetra::Vector<> >& vec )
   : vec_( vec )
   {}

   void setValue( uint_t idx, real_t value ) { vec_->replaceGlobalValue( idx, value ); }

   real_t getValue( uint_t idx ) const
   {
      const auto             map        = vec_->getMap();
      const auto             localIndex = map->getLocalElement( idx );
      return vec_->getData().operator[]( localIndex );
   }

 private:
   RCP< Tpetra::Vector<> > vec_;
};

} // namespace hyteg