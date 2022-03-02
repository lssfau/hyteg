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

#include "hyteg/composites/P1P0StokesFunction.hpp"
#include "hyteg/mixedoperators/P0ScalarToP1VectorOperator.hpp"
#include "hyteg/mixedoperators/P1VectorToP0ScalarOperator.hpp"
#include "hyteg/mixedoperators/P1ToP0Operator.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"

namespace hyteg {

class P1P0StokesOperator : public Operator< P1P0StokesFunction< real_t >, P1P0StokesFunction< real_t > >
{
 public:
   typedef P1ConstantVectorLaplaceOperator VelocityBlockOperator_T;
   typedef P1ConstantLaplaceOperator       VelocityOperator_T;
   // typedef P1P0StokesBlockPreconditioner BlockPreconditioner_T;

   P1P0StokesOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , Lapl( storage, minLevel, maxLevel )
   // , div( storage, minLevel, maxLevel )
   , divT( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   {}

   void apply( const P1P0StokesFunction< real_t >& src,
               const P1P0StokesFunction< real_t >& dst,
               const uint_t                        level,
               const DoFType                       flag ) const
   {
      Lapl.apply( src.uvw(), dst.uvw(), level, flag, Replace );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      // div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1P0StokesFunction< idx_t >&          src,
                  const P1P0StokesFunction< idx_t >&          dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      Lapl.toMatrix( mat, src.uvw(), dst.uvw(), level, flag );
      divT.toMatrix( mat, src.p(), dst.uvw(), level, flag );
      // div.toMatrix( mat, src.uvw(), dst.p(), level, flag );
   }

   VelocityBlockOperator_T    Lapl;
   // P1ToP0ConstantDivOperator  div;
   P0ToP1ConstantDivTOperator divT;

   bool hasGlobalCells_;
};

} // namespace hyteg
