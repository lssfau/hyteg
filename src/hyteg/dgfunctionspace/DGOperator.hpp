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

#include "hyteg/celldofspace/CellDoFIndexing.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/solvers/Smoothables.hpp"

namespace hyteg {
namespace dg {

using walberla::int_c;
using walberla::real_t;

class DGOperator : public Operator< DGFunction< real_t >, DGFunction< real_t > >,
                   public WeightedJacobiSmoothable< DGFunction< real_t > >,
                   public OperatorWithInverseDiagonal< DGFunction< real_t > >
{
 public:
   DGOperator( const std::shared_ptr< PrimitiveStorage >& storage,
               uint_t                                     minLevel,
               uint_t                                     maxLevel,
               const std::shared_ptr< DGForm >&           form );

   void apply( const DGFunction< real_t >& src,
               const DGFunction< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType ) const override;

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >&             mat,
                  const typename srcType::template FunctionType< idx_t >& src,
                  const typename dstType::template FunctionType< idx_t >& dst,
                  size_t                                                  level,
                  DoFType                                                 flag ) const override;

   void smooth_jac( const DGFunction< real_t >& dst,
                    const DGFunction< real_t >& rhs,
                    const DGFunction< real_t >& tmp,
                    real_t                      relax,
                    size_t                      level,
                    DoFType                     flag ) const override;

   [[nodiscard]] std::shared_ptr< DGFunction< real_t > > getInverseDiagonalValues() const override;

 private:
   void applyInner( const DGFunction< real_t >& src,
                    const DGFunction< real_t >& dst,
                    size_t                      level,
                    DoFType                     flag,
                    UpdateType                  updateType = Replace ) const;

   void toMatrixInner( const std::shared_ptr< SparseMatrixProxy >&             mat,
                       const typename srcType::template FunctionType< idx_t >& src,
                       const typename dstType::template FunctionType< idx_t >& dst,
                       size_t                                                  level,
                       DoFType                                                 flag ) const;

   std::shared_ptr< DGForm > form_;
};

} // namespace dg
} // namespace hyteg