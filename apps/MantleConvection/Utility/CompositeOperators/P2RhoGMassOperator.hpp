/*
 * Copyright (c) 2023-2025 Andreas Burkhart.
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
#include "../../Operators/rho_g_mass/P2ElementwiseRhoGMassAnnulusMap_0_0.hpp"
#include "../../Operators/rho_g_mass/P2ElementwiseRhoGMassAnnulusMap_1_0.hpp"
#include "../../Operators/rho_g_mass/P2ElementwiseRhoGMassIcosahedralShellMap_0_0.hpp"
#include "../../Operators/rho_g_mass/P2ElementwiseRhoGMassIcosahedralShellMap_1_0.hpp"
#include "../../Operators/rho_g_mass/P2ElementwiseRhoGMassIcosahedralShellMap_2_0.hpp"
#include "../../Operators/rho_g_mass/P2ElementwiseRhoGMass_0_0.hpp"
#include "../../Operators/rho_g_mass/P2ElementwiseRhoGMass_1_0.hpp"
#include "../../Operators/rho_g_mass/P2ElementwiseRhoGMass_2_0.hpp"

namespace hyteg {

namespace mcoperators {

using walberla::real_t;

class P2RhoGMassOperator : public Operator< P2Function< real_t >, P2VectorFunction< real_t > >
{
 public:
   P2RhoGMassOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                       size_t                                     minLevel,
                       size_t                                     maxLevel,
                       P1Function< real_t >&                      rho )
   : Operator( storage, minLevel, maxLevel )
   , storage_( storage )
   {
      rhoMassOp0_ = std::make_shared< mcoperators::P2ElementwiseRhoGMass_0_0 >( storage, minLevel, maxLevel, rho );
      rhoMassOp1_ = std::make_shared< mcoperators::P2ElementwiseRhoGMass_1_0 >( storage, minLevel, maxLevel, rho );

      if ( storage_->hasGlobalCells() )
      {
         rhoMassOp2_ = std::make_shared< mcoperators::P2ElementwiseRhoGMass_2_0 >( storage, minLevel, maxLevel, rho );
      }
   }

   void apply( const P2Function< real_t >&       src,
               const P2VectorFunction< real_t >& dst,
               const uint_t                      level,
               const DoFType                     flag,
               const UpdateType                  updateType = Replace ) const
   {
      rhoMassOp0_->apply( src, dst.component( 0 ), level, flag, updateType );
      rhoMassOp1_->apply( src, dst.component( 1 ), level, flag, updateType );

      if ( storage_->hasGlobalCells() )
      {
         rhoMassOp2_->apply( src, dst.component( 2 ), level, flag, updateType );
      }
   }

 private:
   std::shared_ptr< mcoperators::P2ElementwiseRhoGMass_0_0 > rhoMassOp0_;
   std::shared_ptr< mcoperators::P2ElementwiseRhoGMass_1_0 > rhoMassOp1_;
   std::shared_ptr< mcoperators::P2ElementwiseRhoGMass_2_0 > rhoMassOp2_;

   const std::shared_ptr< PrimitiveStorage > storage_;
};

class P2RhoGMassIcosahedralShellMapOperator : public Operator< P2Function< real_t >, P2VectorFunction< real_t > >
{
 public:
   P2RhoGMassIcosahedralShellMapOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                          size_t                                     minLevel,
                                          size_t                                     maxLevel,
                                          P1Function< real_t >&                      rho )
   : Operator( storage, minLevel, maxLevel )
   , storage_( storage )
   {
      rhoMassOp0_ = std::make_shared< mcoperators::P2ElementwiseRhoGMassIcosahedralShellMap_0_0 >(
          storage, minLevel, maxLevel, rho );
      rhoMassOp1_ = std::make_shared< mcoperators::P2ElementwiseRhoGMassIcosahedralShellMap_1_0 >(
          storage, minLevel, maxLevel, rho );
      rhoMassOp2_ = std::make_shared< mcoperators::P2ElementwiseRhoGMassIcosahedralShellMap_2_0 >(
          storage, minLevel, maxLevel, rho );

      if ( !storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "P2RhoGMassAnnulusMapOperator in 2D?" );
      }
   }

   void apply( const P2Function< real_t >&       src,
               const P2VectorFunction< real_t >& dst,
               const uint_t                      level,
               const DoFType                     flag,
               const UpdateType                  updateType = Replace ) const
   {
      rhoMassOp0_->apply( src, dst.component( 0 ), level, flag, updateType );
      rhoMassOp1_->apply( src, dst.component( 1 ), level, flag, updateType );
      rhoMassOp2_->apply( src, dst.component( 2 ), level, flag, updateType );
   }

 private:
   std::shared_ptr< mcoperators::P2ElementwiseRhoGMassIcosahedralShellMap_0_0 > rhoMassOp0_;
   std::shared_ptr< mcoperators::P2ElementwiseRhoGMassIcosahedralShellMap_1_0 > rhoMassOp1_;
   std::shared_ptr< mcoperators::P2ElementwiseRhoGMassIcosahedralShellMap_2_0 > rhoMassOp2_;

   const std::shared_ptr< PrimitiveStorage > storage_;
};

class P2RhoGMassAnnulusMapOperator : public Operator< P2Function< real_t >, P2VectorFunction< real_t > >
{
 public:
   P2RhoGMassAnnulusMapOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                 size_t                                     minLevel,
                                 size_t                                     maxLevel,
                                 P1Function< real_t >&                      rho )
   : Operator( storage, minLevel, maxLevel )
   , storage_( storage )
   {
      rhoMassOp0_ =
          std::make_shared< mcoperators::P2ElementwiseRhoGMassAnnulusMap_0_0 >( storage, minLevel, maxLevel, rho );
      rhoMassOp1_ =
          std::make_shared< mcoperators::P2ElementwiseRhoGMassAnnulusMap_1_0 >( storage, minLevel, maxLevel, rho );

      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "P2RhoGMassAnnulusMapOperator in 3D?" );
      }
   }

   void apply( const P2Function< real_t >&       src,
               const P2VectorFunction< real_t >& dst,
               const uint_t                      level,
               const DoFType                     flag,
               const UpdateType                  updateType = Replace ) const
   {
      rhoMassOp0_->apply( src, dst.component( 0 ), level, flag, updateType );
      rhoMassOp1_->apply( src, dst.component( 1 ), level, flag, updateType );
   }

 private:
   std::shared_ptr< mcoperators::P2ElementwiseRhoGMassAnnulusMap_0_0 > rhoMassOp0_;
   std::shared_ptr< mcoperators::P2ElementwiseRhoGMassAnnulusMap_1_0 > rhoMassOp1_;

   const std::shared_ptr< PrimitiveStorage > storage_;
};

} // namespace mcoperators

} // namespace hyteg
