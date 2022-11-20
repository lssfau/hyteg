/*
* Copyright (c) 2022 Andreas Wagner.
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

#include <utility>

#include "hyteg/celldofspace/CellDoFIndexing.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/indexing/Common.hpp"
#include "hyteg/indexing/MacroCellIndexing.hpp"
#include "hyteg/indexing/MacroEdgeIndexing.hpp"
#include "hyteg/indexing/MacroFaceIndexing.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/mixedoperators/DGToP1Operator.hpp"
#include "hyteg/mixedoperators/P1ToDGOperator.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/solvers/Smoothables.hpp"

#include "DGOperator.hpp"

namespace hyteg {
namespace dg {

using facedof::FaceType;
using indexing::Index;
using volumedofspace::indexing::VolumeDoFMemoryLayout;
using walberla::int_c;
using walberla::real_t;

class P1WithDGFormOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
 public:
   P1WithDGFormOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                         uint_t                                     minLevel,
                         uint_t                                     maxLevel,
                         std::shared_ptr< DGForm >                  form )
   : Operator< P1Function< real_t >, P1Function< real_t > >( storage, minLevel, maxLevel )
   , form_( std::move( form ) )
   , srcIdx_( "src", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )
   , dstIdx_( "dst", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )
   , srcReal_( "src", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )
   , dstReal_( "dst", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )
   , opMain( storage, minLevel, maxLevel, form_ )
   , opP1ToDGIdx( storage, minLevel, maxLevel, std::make_shared< P1ToDG1InterpolationForm >() )
   , opP1ToDGReal( storage, minLevel, maxLevel, std::make_shared< P1ToDG1InterpolationForm >() )
   , opDGToP1Real( storage, minLevel, maxLevel, std::make_shared< P1ToDG1InterpolationForm >() )
   {}

   void apply( const P1Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType ) const override
   {
      opP1ToDGReal.apply( src, srcReal_, level, All, Replace );

      if ( updateType != Replace )
         opP1ToDGReal.apply( dst, dstReal_, level, All, Replace );
      opMain.apply( srcReal_, dstReal_, level, All, updateType );

      opDGToP1Real.apply( dstReal_, dst, level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1Function< idx_t >&                  src,
                  const P1Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      opP1ToDGIdx.apply( src, srcIdx_, level, All, Replace );
      opP1ToDGIdx.apply( dst, dstIdx_, level, All, Replace );
      opMain.toMatrix( mat, srcIdx_, dstIdx_, level, flag );
   }

 private:
   std::shared_ptr< DGForm > form_;

   DGFunction< idx_t > srcIdx_;
   DGFunction< idx_t > dstIdx_;

   DGFunction< real_t > srcReal_;
   DGFunction< real_t > dstReal_;

   DGOperator opMain;

   P1ToDGOperator< P1ToDG1InterpolationForm, idx_t >  opP1ToDGIdx;
   P1ToDGOperator< P1ToDG1InterpolationForm, real_t > opP1ToDGReal;
   DGToP1Operator< P1ToDG1InterpolationForm, real_t > opDGToP1Real;
};

} // namespace dg
} // namespace hyteg