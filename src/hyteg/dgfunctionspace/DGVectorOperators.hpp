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

#include "hyteg/dgfunctionspace/DGVectorLaplaceForm.hpp"
#include "hyteg/dgfunctionspace/DGVectorMassForm.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/dgfunctionspace/DGVectorFunction.hpp"
#include "hyteg/operators/VectorToVectorOperator.hpp"

namespace hyteg {
namespace dg {

template < typename ValueType >
class DGVectorLaplaceOperator : public VectorToVectorOperator< ValueType, DGVectorFunction, DGVectorFunction >
{
 public:
   DGVectorLaplaceOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : VectorToVectorOperator< ValueType, DGVectorFunction, DGVectorFunction >( storage, minLevel, maxLevel )
   {
      if ( this->dim_ == 3 )
      {
         WALBERLA_ABORT("not implemented yet.");
      }
      else
      {
         this->subOper_[0][0] = std::make_shared< DGOperator > ( storage, minLevel, maxLevel, std::make_shared< forms::DGVectorLaplaceFormP1P1_00 >() );
         this->subOper_[0][1] = std::make_shared< DGOperator > ( storage, minLevel, maxLevel, std::make_shared< forms::DGVectorLaplaceFormP1P1_01 >() );
         this->subOper_[1][0] = std::make_shared< DGOperator > ( storage, minLevel, maxLevel, std::make_shared< forms::DGVectorLaplaceFormP1P1_10 >() );
         this->subOper_[1][1] = std::make_shared< DGOperator > ( storage, minLevel, maxLevel, std::make_shared< forms::DGVectorLaplaceFormP1P1_11 >() );
      }
   }
};

template < typename ValueType >
class DGVectorMassOperator : public VectorToVectorOperator< ValueType, DGVectorFunction, DGVectorFunction >
{
 public:
   DGVectorMassOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
       : VectorToVectorOperator< ValueType, DGVectorFunction, DGVectorFunction >( storage, minLevel, maxLevel )
   {
      if ( this->dim_ == 3 )
      {
         WALBERLA_ABORT("not implemented yet.");
      }
      else
      {
         this->subOper_[0][0] = std::make_shared< DGOperator > ( storage, minLevel, maxLevel, std::make_shared< forms::DGVectorMassFormP1P1_00 >() );
         this->subOper_[0][1] = std::make_shared< DGOperator > ( storage, minLevel, maxLevel, std::make_shared< forms::DGVectorMassFormP1P1_01 >() );
         this->subOper_[1][0] = std::make_shared< DGOperator > ( storage, minLevel, maxLevel, std::make_shared< forms::DGVectorMassFormP1P1_10 >() );
         this->subOper_[1][1] = std::make_shared< DGOperator > ( storage, minLevel, maxLevel, std::make_shared< forms::DGVectorMassFormP1P1_11 >() );
      }
   }
};

}
}