/*
 * Copyright (c) 2024 Ponsuganth Ilangovan P.
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

// #include "hyteg/operatorgeneration/generated/EpsilonRotation/P2VectorElementwiseEpsilonRotationAbsIcosahedralShellMapOperator.hpp"
#include "hyteg/p2functionspace/P2RotationOperator.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockPreconditioners.hpp"
#include "hyteg_operators/operators/epsilon/P2VectorElementwiseEpsilonRotationP0ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/divergence/P2VectorToP1ElementwiseDivergenceRotationIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/gradient/P1ToP2VectorElementwiseGradientRotationIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg/operatorgeneration/EpsilonRotation/P2VectorElementwiseEpsilonRotation_IcosahedralShellMap_fused_quadloops_float64.hpp"
#include "hyteg_operators/operators/epsilon/P2VectorElementwiseEpsilonRotationP1ViscosityIcosahedralShellMap.hpp"
#include "hyteg/operatorgeneration/generated/FullStokesRotationWithFSPenalty/P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator.hpp"

namespace hyteg {

class P2ViscousIcosahedralShellMapOperatorRotation : public Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >,
                                                     public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >
{
 public:
   // typedef operatorgeneration::P2VectorElementwiseEpsilonRotationP0ViscosityIcosahedralShellMap ViscousOperator_T;
   // typedef operatorgeneration::P2VectorElementwiseEpsilonRotationP1ViscosityIcosahedralShellMap ViscousOperator_T;
   typedef operatorgeneration::P2VectorElementwiseEpsilonRotation_IcosahedralShellMap_fused_quadloops_float64 ViscousOperator_T;
   // typedef operatorgeneration::P2VectorElementwiseFullStokesRotationWithFusedFSPenaltyVectorizedOperator ViscousOperator_T;
   // using P2VectorElementwiseEpsilonRotationAbsOperator =
   //     operatorgeneration::P2VectorElementwiseEpsilonRotationAbsIcosahedralShellMapOperator;

   P2ViscousIcosahedralShellMapOperatorRotation( const std::shared_ptr< PrimitiveStorage >& storage,
                                                 uint_t                                     minLevel,
                                                 uint_t                                     maxLevel,
                                                 P0Function< real_t >&                      mu,
                                                 P2Function< real_t >&                      nx,
                                                 P2Function< real_t >&                      ny,
                                                 P2Function< real_t >&                      nz,
                                                 real_t                                     rotFactor )
   : Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >( storage, minLevel, maxLevel )
   , tmp_( "tmp__P2ViscousIcosahedralShellMapOperatorRotation", storage, minLevel, maxLevel )
   , ones_( "ones__P2ViscousIcosahedralShellMapOperatorRotation", storage, minLevel, maxLevel )
   , viscousOperator( storage, minLevel, maxLevel, mu, nx, ny, nz, rotFactor )
   // , epsilonAbsOperator( storage, minLevel, maxLevel, mu, nx, ny, nz )
   {
      invDiag_ = std::make_shared< P2VectorFunction< real_t > >(
          "invDiag__P2ViscousIcosahedralShellMapOperatorRotation", storage, minLevel, maxLevel );
   }

   void apply( const P2VectorFunction< real_t >& src,
               const P2VectorFunction< real_t >& dst,
               const uint_t                      level,
               const DoFType                     flag,
               const UpdateType                  updateType = Replace ) const
   {
      viscousOperator.apply( src, dst, level, flag, updateType );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2VectorFunction< idx_t >&            numeratorSrc,
                  const P2VectorFunction< idx_t >&            numeratorDst,
                  uint_t                                      level,
                  DoFType                                     flag ) const
   {
      viscousOperator.toMatrix( mat, numeratorSrc, numeratorDst, level, flag );
   }

   void computeInverseDiagonalOperatorValues()
   {
      viscousOperator.computeInverseDiagonalOperatorValues();
      // auto viscousInvDiag = viscousOperator.getInverseDiagonalValues();

      // std::function< real_t( const Point3D&, const std::vector< real_t >& ) > invHelper =
      //     []( const Point3D&, const std::vector< real_t >& vals ) { return 1.0 / vals[0]; };

      // for( uint_t level = tmp_.getMinLevel(); level <= tmp_.getMaxLevel(); level++ )
      // {
      //    tmp_.component(0u).interpolate( invHelper, { viscousInvDiag->component(0u) }, level, All );
      //    tmp_.component(1u).interpolate( invHelper, { viscousInvDiag->component(1u) }, level, All );
      //    tmp_.component(2u).interpolate( invHelper, { viscousInvDiag->component(2u) }, level, All );

      //    ones_.interpolate( 1.0, level, All );

      //    epsilonAbsOperator.apply( ones_, tmp_, level, All, Add );

      //    invDiag_->component(0u).interpolate( invHelper, { tmp_.component(0u) }, level, All );
      //    invDiag_->component(1u).interpolate( invHelper, { tmp_.component(1u) }, level, All );
      //    invDiag_->component(2u).interpolate( invHelper, { tmp_.component(2u) }, level, All );
      // }
   }

   std::shared_ptr< P2VectorFunction< real_t > > getInverseDiagonalValues() const
   {
      return viscousOperator.getInverseDiagonalValues();
      // return invDiag_;
   }

   P2VectorFunction< real_t > tmp_;
   P2VectorFunction< real_t > ones_;

   std::shared_ptr< P2VectorFunction< real_t > > invDiag_;

   ViscousOperator_T viscousOperator;

   // P2VectorElementwiseEpsilonRotationAbsOperator epsilonAbsOperator;
};

class P2P1StokesOpgenRotationWrapper : public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   // typedef operatorgeneration::P2VectorElementwiseEpsilonRotationIcosahedralShellMap        StokesViscousOperatorRotationOpgen;
   typedef P2ViscousIcosahedralShellMapOperatorRotation                                     StokesViscousOperatorRotationOpgen;
   typedef operatorgeneration::P1ToP2VectorElementwiseGradientRotationIcosahedralShellMap   StokesGradientOperator;
   typedef operatorgeneration::P2VectorToP1ElementwiseDivergenceRotationIcosahedralShellMap StokesDivergenceOperator;
   typedef StokesViscousOperatorRotationOpgen                                               VelocityOperator_T;

   using StabilizationOperator_T = ZeroOperator< P1Function< real_t >, P1Function< real_t > >;
   using DivOperator_T           = StokesDivergenceOperator;
   using GradOperator_T          = StokesGradientOperator;
   using StabOperator_T          = StabilizationOperator_T;

   P2P1StokesOpgenRotationWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                                   uint_t                                     minLevel,
                                   uint_t                                     maxLevel,
                                   P0Function< real_t >&                      mu,
                                   //   P1Function< real_t >&                      muInv,
                                   P2Function< real_t >& nx,
                                   P2Function< real_t >& ny,
                                   P2Function< real_t >& nz,
                                   real_t                rotFactor,
                                   P2RotationOperator&   rotationOperator,
                                   BoundaryCondition     bcVelocity )
   : Operator( storage, minLevel, maxLevel )
   , stokesViscousOperator_( storage, minLevel, maxLevel, mu, nx, ny, nz, rotFactor )
   , rotationOperator_( rotationOperator )
   , divT( storage, minLevel, maxLevel, nx, ny, nz )
   , div( storage, minLevel, maxLevel, nx, ny, nz )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , stabOp_( storage, minLevel, maxLevel )
   , tmp_( "tmp__P2P1StokesRotationWrapper", storage, minLevel, maxLevel, bcVelocity )
   , tmpdst_( "tmpdst__P2P1StokesRotationWrapper", storage, minLevel, maxLevel, bcVelocity )
   , tmpAssembly_( "tmpAssembly__P2P1StokesRotationWrapper", storage, minLevel, maxLevel, bcVelocity )
   {
      tmpAssembly_.enumerate( maxLevel );
   }

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               uint_t                                  level,
               DoFType                                 flag ) const
   {
      stokesViscousOperator_.apply( src.uvw(), dst.uvw(), level, flag );
      divT.apply( src.p(), dst.uvw(), level, flag, Add );
      div.apply( src.uvw(), dst.p(), level, flag, Replace );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< idx_t >&      numeratorSrc,
                  const P2P1TaylorHoodFunction< idx_t >&      numeratorDst,
                  uint_t                                      level,
                  DoFType                                     flag ) const
   {
      stokesViscousOperator_.viscousOperator.toMatrix( mat, numeratorSrc.uvw(), numeratorDst.uvw(), level, flag );
      divT.toMatrix( mat, numeratorSrc.p(), numeratorDst.uvw(), level, flag );
      div.toMatrix( mat, numeratorSrc.uvw(), numeratorDst.p(), level, flag );
   }

   const VelocityOperator_T& getA() const { return stokesViscousOperator_; }
   const DivOperator_T&      getB() const { return div; }
   const GradOperator_T&     getBT() const { return divT; }
   const StabOperator_T&     getStab() const { return stabOp_; }

   VelocityOperator_T& getA() { return stokesViscousOperator_; }
   DivOperator_T&      getB() { return div; }
   GradOperator_T&     getBT() { return divT; }
   StabOperator_T&     getStab() { return stabOp_; }

   StokesViscousOperatorRotationOpgen stokesViscousOperator_;
   P2RotationOperator&                rotationOperator_;

   StokesGradientOperator   divT;
   StokesDivergenceOperator div;

   // operatorgeneration::P1ElementwiseKMassIcosahedralShellMap pspg_inv_diag_;
   P1PSPGInvDiagOperator pspg_inv_diag_;

   StabilizationOperator_T stabOp_;

   P2P1TaylorHoodFunction< real_t > tmp_;
   P2P1TaylorHoodFunction< real_t > tmpdst_;

   P2P1TaylorHoodFunction< idx_t > tmpAssembly_;
};

} // namespace hyteg