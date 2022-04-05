/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1dgefunctionspace/P1DGEFunction.hpp"
#include "hyteg/dgfunctionspace/DGVectorMassForm.hpp"

namespace hyteg {

template < typename ValueType >
class P1DGEMassOperator final : public Operator< P1DGEFunction< real_t >, P1DGEFunction< real_t > > {
 public:
   P1DGEMassOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                  uint_t                                     minLevel,
                  uint_t                                     maxLevel);

   void apply( const P1DGEFunction< real_t >& src,
               const P1DGEFunction< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType ) const override
   {

   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1DGEFunction< idx_t >&                  src,
                  const P1DGEFunction< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {

   }

 private:


};

}
