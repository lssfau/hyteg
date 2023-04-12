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

#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/egfunctionspace/EGEpsilonFormNitscheBC.hpp"
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

class ScalarP1WithDGFormOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
 public:
   ScalarP1WithDGFormOperator( const std::shared_ptr< PrimitiveStorage >& storage,
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
      //srcReal_.interpolate( 0, level, All );
      opP1ToDGReal.apply( src, srcReal_, level, flag, Replace );

      //if ( updateType != Replace )
      //   opP1ToDGReal.apply( dst, dstReal_, level, flag, Replace );

      //dstReal_.interpolate( 0, level, All );
      opMain.apply( srcReal_, dstReal_, level, flag, Replace );

      opDGToP1Real.apply( dstReal_, dst, level, flag, updateType );
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

class VectorialP1WithEpsilonFormOperator : public Operator< P1VectorFunction< real_t >, P1VectorFunction< real_t > >
{
 public:
   VectorialP1WithEpsilonFormOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                       uint_t                                     minLevel,
                                       uint_t                                     maxLevel,
                                       std::function< real_t( const Point3D& ) >  viscosity )
   : Operator< P1VectorFunction< real_t >, P1VectorFunction< real_t > >( storage, minLevel, maxLevel )
   , srcIdx_0( "src_0", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )
   , srcIdx_1( "src_1", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )
   , srcIdx_2( "src_2", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )

   , dstIdx_0( "dst_0", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )
   , dstIdx_1( "dst_1", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )
   , dstIdx_2( "dst_2", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )

   , srcReal_0( "src_0", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )
   , srcReal_1( "src_1", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )
   , srcReal_2( "src_2", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )

   , dstReal_0( "dst_0", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )
   , dstReal_1( "dst_1", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )
   , dstReal_2( "dst_2", storage, minLevel, maxLevel, std::make_shared< DGBasisLinearLagrange_Example >(), 1 )

   , eps_00( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_00 >( viscosity ) )
   , eps_01( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_01 >( viscosity ) )
   , eps_02( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_02 >( viscosity ) )
   , eps_10( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_10 >( viscosity ) )
   , eps_11( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_11 >( viscosity ) )
   , eps_12( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_12 >( viscosity ) )
   , eps_20( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_20 >( viscosity ) )
   , eps_21( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_21 >( viscosity ) )
   , eps_22( storage, minLevel, maxLevel, std::make_shared< eg::EGEpsilonFormNitscheBC_P1P1_22 >( viscosity ) )

   , opP1ToDGIdx( storage, minLevel, maxLevel, std::make_shared< P1ToDG1InterpolationForm >() )
   , opP1ToDGReal( storage, minLevel, maxLevel, std::make_shared< P1ToDG1InterpolationForm >() )
   , opDGToP1Real( storage, minLevel, maxLevel, std::make_shared< P1ToDG1InterpolationForm >() )
   {
      if ( !storage->hasGlobalCells() )
         WALBERLA_ABORT( "Not implemented." )
   }

   void apply( const P1VectorFunction< real_t >& src,
               const P1VectorFunction< real_t >& dst,
               size_t                            level,
               DoFType                           flag,
               UpdateType                        updateType ) const override
   {
      // interpolate p1 to dg function
      opP1ToDGReal.apply( src.component( 0 ), srcReal_0, level, All, Replace );
      opP1ToDGReal.apply( src.component( 1 ), srcReal_1, level, All, Replace );
      opP1ToDGReal.apply( src.component( 2 ), srcReal_2, level, All, Replace );

      // transfer dst to dg if we e.g. add to it
      if ( updateType != Replace )
      {
         opP1ToDGReal.apply( dst.component( 0 ), dstReal_0, level, All, Replace );
         opP1ToDGReal.apply( dst.component( 1 ), dstReal_1, level, All, Replace );
         opP1ToDGReal.apply( dst.component( 2 ), dstReal_2, level, All, Replace );
      }

      // apply dg operator with weak bcs
      eps_00.apply( srcReal_0, dstReal_0, level, All, updateType );
      eps_01.apply( srcReal_1, dstReal_0, level, All, Add );
      eps_02.apply( srcReal_2, dstReal_0, level, All, Add );

      eps_10.apply( srcReal_0, dstReal_1, level, All, updateType );
      eps_11.apply( srcReal_1, dstReal_1, level, All, Add );
      eps_12.apply( srcReal_2, dstReal_1, level, All, Add );

      eps_20.apply( srcReal_0, dstReal_2, level, All, updateType );
      eps_21.apply( srcReal_1, dstReal_2, level, All, Add );
      eps_22.apply( srcReal_2, dstReal_2, level, All, Add );

      // interpolate back from dg to p1
      opDGToP1Real.apply( dstReal_0, dst.component( 0 ), level, flag, Replace );
      opDGToP1Real.apply( dstReal_1, dst.component( 1 ), level, flag, Replace );
      opDGToP1Real.apply( dstReal_2, dst.component( 2 ), level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1VectorFunction< idx_t >&            src,
                  const P1VectorFunction< idx_t >&            dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      opP1ToDGIdx.apply( src.component( 0 ), srcIdx_0, level, All, Replace );
      opP1ToDGIdx.apply( src.component( 1 ), srcIdx_1, level, All, Replace );
      opP1ToDGIdx.apply( src.component( 2 ), srcIdx_2, level, All, Replace );

      opP1ToDGIdx.apply( dst.component( 0 ), dstIdx_0, level, All, Replace );
      opP1ToDGIdx.apply( dst.component( 1 ), dstIdx_1, level, All, Replace );
      opP1ToDGIdx.apply( dst.component( 2 ), dstIdx_2, level, All, Replace );

      eps_00.toMatrix( mat, srcIdx_0, dstIdx_0, level, flag );
      eps_01.toMatrix( mat, srcIdx_1, dstIdx_0, level, flag );
      eps_02.toMatrix( mat, srcIdx_2, dstIdx_0, level, flag );
      eps_10.toMatrix( mat, srcIdx_0, dstIdx_1, level, flag );
      eps_11.toMatrix( mat, srcIdx_1, dstIdx_1, level, flag );
      eps_12.toMatrix( mat, srcIdx_2, dstIdx_1, level, flag );
      eps_20.toMatrix( mat, srcIdx_0, dstIdx_2, level, flag );
      eps_21.toMatrix( mat, srcIdx_1, dstIdx_2, level, flag );
      eps_22.toMatrix( mat, srcIdx_2, dstIdx_2, level, flag );
   }

 private:
   DGFunction< idx_t > srcIdx_0;
   DGFunction< idx_t > srcIdx_1;
   DGFunction< idx_t > srcIdx_2;

   DGFunction< idx_t > dstIdx_0;
   DGFunction< idx_t > dstIdx_1;
   DGFunction< idx_t > dstIdx_2;

   DGFunction< real_t > srcReal_0;
   DGFunction< real_t > srcReal_1;
   DGFunction< real_t > srcReal_2;

   DGFunction< real_t > dstReal_0;
   DGFunction< real_t > dstReal_1;
   DGFunction< real_t > dstReal_2;

   DGOperator eps_00;
   DGOperator eps_01;
   DGOperator eps_02;
   DGOperator eps_10;
   DGOperator eps_11;
   DGOperator eps_12;
   DGOperator eps_20;
   DGOperator eps_21;
   DGOperator eps_22;

   P1ToDGOperator< P1ToDG1InterpolationForm, idx_t >  opP1ToDGIdx;
   P1ToDGOperator< P1ToDG1InterpolationForm, real_t > opP1ToDGReal;
   DGToP1Operator< P1ToDG1InterpolationForm, real_t > opDGToP1Real;
};
} // namespace dg
} // namespace hyteg