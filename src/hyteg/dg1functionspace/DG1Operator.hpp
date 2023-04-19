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

#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/solvers/Smoothables.hpp"

namespace hyteg {

using walberla::int_c;
using walberla::real_t;

template < typename DG1Form >
class DG1Operator : public Operator< DG1Function< real_t >, DG1Function< real_t > >
{
 public:
   DG1Operator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : DG1Operator( storage, minLevel, maxLevel, std::make_shared< DG1Form >() )
   {}

   DG1Operator( const std::shared_ptr< PrimitiveStorage >& storage,
                uint_t                                     minLevel,
                uint_t                                     maxLevel,
                const std::shared_ptr< DG1Form >&          form )
   : Operator< DG1Function< real_t >, DG1Function< real_t > >( storage, minLevel, maxLevel )
   , dgOperator_( storage, minLevel, maxLevel, form )
   {}

   void apply( const DG1Function< real_t >& src,
               const DG1Function< real_t >& dst,
               size_t                       level,
               DoFType                      flag,
               UpdateType                   updateType ) const override
   {
      dgOperator_.apply( *src.getDGFunction(), *dst.getDGFunction(), level, flag, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const DG1Function< idx_t >&                 src,
                  const DG1Function< idx_t >&                 dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      dgOperator_.toMatrix( mat, *src.getDGFunction(), *dst.getDGFunction(), level, flag );
   }

 private:
   dg::DGOperator dgOperator_;
};

} // namespace hyteg