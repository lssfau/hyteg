/*
 * Copyright (c) 2025 Andreas Burkhart.
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
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergenceAnnulusMap_0_0.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergenceAnnulusMap_0_1.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_0.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_1.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_2.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergence_0_0.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergence_0_1.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergence_0_2.hpp"

namespace hyteg {

namespace mcoperators {

using walberla::real_t;

class P2ToP1GradRhoRhoDivergenceOperator : public Operator< P2VectorFunction< real_t >, P1Function< real_t > >
{
 public:
   P2ToP1GradRhoRhoDivergenceOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                             size_t                                     minLevel,
                             size_t                                     maxLevel,
                             P1Function< real_t >&                      invRho,
                             P1Function< real_t >&                      rho )
   : Operator( storage, minLevel, maxLevel )
   , storage_( storage )
   {
      GradRhoRhoOp0_ =
          std::make_shared< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergence_0_0 >( storage, minLevel, maxLevel, invRho, rho );
      GradRhoRhoOp1_ =
          std::make_shared< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergence_0_1 >( storage, minLevel, maxLevel, invRho, rho );
      if ( storage_->hasGlobalCells() )
      {
         GradRhoRhoOp2_ =
             std::make_shared< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergence_0_2 >( storage, minLevel, maxLevel, invRho, rho );
      }
   }

   void apply( const P2VectorFunction< real_t >& src,
               const P1Function< real_t >&       dst,
               const uint_t                      level,
               const DoFType                     flag,
               const UpdateType                  updateType = Replace ) const
   {
      GradRhoRhoOp0_->apply( src.component( 0 ), dst, level, flag, updateType );
      GradRhoRhoOp1_->apply( src.component( 1 ), dst, level, flag, Add );

      if ( storage_->hasGlobalCells() )
      {
         GradRhoRhoOp2_->apply( src.component( 2 ), dst, level, flag, Add );
      }
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2VectorFunction< idx_t >&            src,
                  const P1Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      GradRhoRhoOp0_->toMatrix( mat, src[0], dst, level, flag );
      GradRhoRhoOp1_->toMatrix( mat, src[1], dst, level, flag );
      if ( storage_->hasGlobalCells() )
      {
         GradRhoRhoOp2_->toMatrix( mat, src[2], dst, level, flag );
      }
   }

 private:
   std::shared_ptr< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergence_0_0 > GradRhoRhoOp0_;
   std::shared_ptr< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergence_0_1 > GradRhoRhoOp1_;
   std::shared_ptr< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergence_0_2 > GradRhoRhoOp2_;

   const std::shared_ptr< PrimitiveStorage > storage_;
};

class P2ToP1GradRhoRhoDivergenceIcosahedralShellMapOperator : public Operator< P2VectorFunction< real_t >, P1Function< real_t > >
{
 public:
   P2ToP1GradRhoRhoDivergenceIcosahedralShellMapOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                size_t                                     minLevel,
                                                size_t                                     maxLevel,
                                                P1Function< real_t >&                      invRho,
                                                P1Function< real_t >&                      rho )
   : Operator( storage, minLevel, maxLevel )
   , storage_( storage )
   {
      GradRhoRhoOp0_ = std::make_shared< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_0 >(
          storage, minLevel, maxLevel, invRho, rho );
      GradRhoRhoOp1_ = std::make_shared< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_1 >(
          storage, minLevel, maxLevel, invRho, rho );
      GradRhoRhoOp2_ = std::make_shared< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_2 >(
          storage, minLevel, maxLevel, invRho, rho );

      if ( !storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "P2ToP1GradRhoRhoDivergenceIcosahedralShellMapOperator in 2D?" );
      }
   }

   void apply( const P2VectorFunction< real_t >& src,
               const P1Function< real_t >&       dst,
               const uint_t                      level,
               const DoFType                     flag,
               const UpdateType                  updateType = Replace ) const
   {
      GradRhoRhoOp0_->apply( src.component( 0 ), dst, level, flag, updateType );
      GradRhoRhoOp1_->apply( src.component( 1 ), dst, level, flag, Add );
      GradRhoRhoOp2_->apply( src.component( 2 ), dst, level, flag, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2VectorFunction< idx_t >&            src,
                  const P1Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      GradRhoRhoOp0_->toMatrix( mat, src[0], dst, level, flag );
      GradRhoRhoOp1_->toMatrix( mat, src[1], dst, level, flag );
      GradRhoRhoOp2_->toMatrix( mat, src[2], dst, level, flag );
   }

 private:
   std::shared_ptr< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_0 > GradRhoRhoOp0_;
   std::shared_ptr< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_1 > GradRhoRhoOp1_;
   std::shared_ptr< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_2 > GradRhoRhoOp2_;

   const std::shared_ptr< PrimitiveStorage > storage_;
};

class P2ToP1GradRhoRhoDivergenceAnnulusMapOperator : public Operator< P2VectorFunction< real_t >, P1Function< real_t > >
{
 public:
   P2ToP1GradRhoRhoDivergenceAnnulusMapOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                       size_t                                     minLevel,
                                       size_t                                     maxLevel,
                                       P1Function< real_t >&                      invRho,
                                       P1Function< real_t >&                      rho )
   : Operator( storage, minLevel, maxLevel )
   , storage_( storage )
   {
      GradRhoRhoOp0_ =
          std::make_shared< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceAnnulusMap_0_0 >( storage, minLevel, maxLevel, invRho, rho );
      GradRhoRhoOp1_ =
          std::make_shared< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceAnnulusMap_0_1 >( storage, minLevel, maxLevel, invRho, rho );

      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "P2ToP1GradRhoRhoDivergenceAnnulusMapOperator in 3D?" );
      }
   }

   void apply( const P2VectorFunction< real_t >& src,
               const P1Function< real_t >&       dst,
               const uint_t                      level,
               const DoFType                     flag,
               const UpdateType                  updateType = Replace ) const
   {
      GradRhoRhoOp0_->apply( src.component( 0 ), dst, level, flag, updateType );
      GradRhoRhoOp1_->apply( src.component( 1 ), dst, level, flag, Add );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2VectorFunction< idx_t >&            src,
                  const P1Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      GradRhoRhoOp0_->toMatrix( mat, src[0], dst, level, flag );
      GradRhoRhoOp1_->toMatrix( mat, src[1], dst, level, flag );
   }

 private:
   std::shared_ptr< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceAnnulusMap_0_0 > GradRhoRhoOp0_;
   std::shared_ptr< mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceAnnulusMap_0_1 > GradRhoRhoOp1_;

   const std::shared_ptr< PrimitiveStorage > storage_;
};

} // namespace mcoperators

} // namespace hyteg