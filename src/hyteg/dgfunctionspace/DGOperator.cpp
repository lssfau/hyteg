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

#include "hyteg/dgfunctionspace/DGOperator.hpp"

#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"

namespace hyteg {
namespace dg {

DGOperator::DGOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                        uint_t                                     minLevel,
                        uint_t                                     maxLevel,
                        const std::shared_ptr< DGForm >&           form )
: Operator< DGFunction< real_t >, DGFunction< real_t > >( storage, minLevel, maxLevel )
, form_( form )
{}

void DGOperator::apply( const DGFunction< real_t >& src,
                        const DGFunction< real_t >& dst,
                        size_t                      level,
                        DoFType                     flag,
                        UpdateType                  updateType ) const
{
   // TODO: communicate

   assembleAndOrApply< real_t >( src, dst, level, flag, nullptr, updateType );
};

void DGOperator::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                           const DGFunction< idx_t >&                  src,
                           const DGFunction< idx_t >&                  dst,
                           size_t                                      level,
                           DoFType                                     flag ) const
{
   // TODO: communicate

   assembleAndOrApply< idx_t >( src, dst, level, flag, mat, Replace );
}

} // namespace dg
} // namespace hyteg
