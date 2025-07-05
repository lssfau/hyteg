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

#include "hyteg/p2functionspace/P2RotationOperator.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockPreconditioners.hpp"
#include "hyteg_operators/operators/divergence/P2VectorToP1ElementwiseDivergenceRotationIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/epsilon/P2VectorElementwiseEpsilonRotationP1ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/epsilon/P2VectorElementwiseEpsilonP1ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/gradient/P1ToP2VectorElementwiseGradientRotationIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"

namespace hyteg {

class P2ViscousIcosahedralShellMapOperatorRotation : public Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >,
                                                     public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >
{
 public:
   typedef operatorgeneration::P2VectorElementwiseEpsilonP1ViscosityIcosahedralShellMap ViscousBaseOperator_T;
   typedef operatorgeneration::P2VectorElementwiseEpsilonRotationP1ViscosityIcosahedralShellMap ViscousOperator_T;

   P2ViscousIcosahedralShellMapOperatorRotation( const std::shared_ptr< PrimitiveStorage >& storage,
                                                 uint_t                                     minLevel,
                                                 uint_t                                     maxLevel,
                                                 const P1Function< real_t >&                mu,
                                                 P2Function< real_t >&                      nx,
                                                 P2Function< real_t >&                      ny,
                                                 P2Function< real_t >&                      nz,
                                                 P2RotationOperator& rotationOperator,
                                                 BoundaryCondition bcVelocity,
                                                 real_t                                     rotFactor )
   : Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >( storage, minLevel, maxLevel )
   , rotationOperator_( rotationOperator )
   , tmp_( "tmp__P2ViscousIcosahedralShellMapOperatorRotation", storage, minLevel, maxLevel, bcVelocity )
   , tmpdst_( "tmpdst__P2ViscousIcosahedralShellMapOperatorRotation", storage, minLevel, maxLevel, bcVelocity )
   , ones_( "ones__P2ViscousIcosahedralShellMapOperatorRotation", storage, minLevel, maxLevel )
   , viscousOperator( storage, minLevel, maxLevel, mu, nx, ny, nz )
   , viscousBaseOperator( storage, minLevel, maxLevel, mu )
   {
      invDiag_ = std::make_shared< P2VectorFunction< real_t > >(
          "invDiag__P2ViscousIcosahedralShellMapOperatorRotation", storage, minLevel, maxLevel );

      viscousOperator.computeInverseDiagonalOperatorValues();
   }

   void apply( const P2VectorFunction< real_t >& src,
               const P2VectorFunction< real_t >& dst,
               const uint_t                      level,
               const DoFType                     flag,
               const UpdateType                  updateType = Replace ) const
   {
      // tmp_.assign( { 1.0 }, { src }, level, All );
      // rotationOperator_.rotate( tmp_, level, FreeslipBoundary, true );

      // {
      //    viscousBaseOperator.apply( tmp_, tmpdst_, level, flag );
      // }

      // rotationOperator_.rotate( tmpdst_, level, FreeslipBoundary );
      // dst.assign( { updateType == Add ? 1.0 : 0.0, 1.0 }, { dst, tmpdst_ }, level, All );

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

   void computeInverseDiagonalOperatorValues() { viscousOperator.computeInverseDiagonalOperatorValues(); }

   std::shared_ptr< P2VectorFunction< real_t > > getInverseDiagonalValues() const
   {
      return viscousOperator.getInverseDiagonalValues();
   }

   P2RotationOperator& rotationOperator_;

   P2VectorFunction< real_t > tmp_;
   P2VectorFunction< real_t > tmpdst_;
   P2VectorFunction< real_t > ones_;

   std::shared_ptr< P2VectorFunction< real_t > > invDiag_;

   ViscousOperator_T viscousOperator;
   ViscousBaseOperator_T viscousBaseOperator;
};

class P2P1StokesOpgenRotationWrapper : public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2ViscousIcosahedralShellMapOperatorRotation                                     StokesViscousOperatorRotationOpgen;
   typedef operatorgeneration::P1ToP2VectorElementwiseGradientRotationIcosahedralShellMap   StokesGradientOperator;
   typedef operatorgeneration::P2VectorToP1ElementwiseDivergenceRotationIcosahedralShellMap StokesDivergenceOperator;
   typedef StokesViscousOperatorRotationOpgen                                               VelocityOperator_T;

   using StabilizationOperator_T = ZeroOperator< P1Function< real_t >, P1Function< real_t > >;
   using DivOperator_T           = StokesDivergenceOperator;
   using GradOperator_T          = StokesGradientOperator;
   using StabOperator_T          = StabilizationOperator_T;

   using StokesBaseOperator_T = operatorgeneration::P2P1StokesEpsilonP1ViscosityIcosahedralShellMapOperator;

   using ViscousOperatorFS_T = StokesViscousOperatorRotationOpgen;
   using SchurOperator_T = operatorgeneration::P1ElementwiseKMassIcosahedralShellMap;

   P2P1StokesOpgenRotationWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                                   uint_t                                     minLevel,
                                   uint_t                                     maxLevel,
                                   const P1Function< real_t >&                mu,
                                   const P1Function< real_t >&                muInv,
                                   P2Function< real_t >&                      nx,
                                   P2Function< real_t >&                      ny,
                                   P2Function< real_t >&                      nz,
                                   real_t                                     rotFactor,
                                   P2RotationOperator&                        rotationOperator,
                                   BoundaryCondition                          bcVelocity )
   : Operator( storage, minLevel, maxLevel )
   , stokesBaseOperator_( storage, minLevel, maxLevel, mu )
   , stokesViscousOperator_( storage, minLevel, maxLevel, mu, nx, ny, nz, rotationOperator, bcVelocity, rotFactor )
   , rotationOperator_( rotationOperator )
   , divT( storage, minLevel, maxLevel, nx, ny, nz )
   , div( storage, minLevel, maxLevel, nx, ny, nz )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , stabOp_( storage, minLevel, maxLevel )
   , schurOperator_( storage, minLevel, maxLevel, muInv )
   , tmp_( "tmp__P2P1StokesRotationWrapper", storage, minLevel, maxLevel, bcVelocity )
   , tmpdst_( "tmpdst__P2P1StokesRotationWrapper", storage, minLevel, maxLevel, bcVelocity )
   , tmpAssembly_( "tmpAssembly__P2P1StokesRotationWrapper", storage, minLevel, maxLevel, bcVelocity )
   {
      tmpAssembly_.enumerate( maxLevel );
      stokesBaseOperator_.getA().computeInverseDiagonalOperatorValues();
      stokesViscousOperator_.computeInverseDiagonalOperatorValues();
      schurOperator_.computeInverseDiagonalOperatorValues();
   }

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               uint_t                                  level,
               DoFType                                 flag,
               const UpdateType                        updateType = Replace ) const override
   {
      tmp_.assign( { 1.0 }, { src }, level, All );
      rotationOperator_.rotate( tmp_, level, FreeslipBoundary, true );
      
      {
         stokesBaseOperator_.getA().apply( tmp_.uvw(), tmpdst_.uvw(), level, flag );
         stokesBaseOperator_.getBT().apply( tmp_.p(), tmpdst_.uvw(), level, flag, Add );
         stokesBaseOperator_.getB().apply( tmp_.uvw(), tmpdst_.p(), level, flag, Replace );
      }

      rotationOperator_.rotate( tmpdst_, level, FreeslipBoundary );
      dst.assign( { 1.0 }, { tmpdst_ }, level, All );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< idx_t >&      numeratorSrc,
                  const P2P1TaylorHoodFunction< idx_t >&      numeratorDst,
                  uint_t                                      level,
                  DoFType                                     flag ) const override
   {
      stokesViscousOperator_.viscousOperator.toMatrix( mat, numeratorSrc.uvw(), numeratorDst.uvw(), level, flag );
      divT.toMatrix( mat, numeratorSrc.p(), numeratorDst.uvw(), level, flag );
      div.toMatrix( mat, numeratorSrc.uvw(), numeratorDst.p(), level, flag );
   }

   const VelocityOperator_T& getA() const { return stokesViscousOperator_; }
   const DivOperator_T&      getB() const { return div; }
   const GradOperator_T&     getBT() const { return divT; }
   const StabOperator_T&     getStab() const { return stabOp_; }
   const SchurOperator_T&    getSchur() const { return schurOperator_; }

   VelocityOperator_T& getA() { return stokesViscousOperator_; }
   DivOperator_T&      getB() { return div; }
   GradOperator_T&     getBT() { return divT; }
   StabOperator_T&     getStab() { return stabOp_; }
   SchurOperator_T&    getSchur() { return schurOperator_; }

   StokesBaseOperator_T stokesBaseOperator_;

   StokesViscousOperatorRotationOpgen stokesViscousOperator_;
   P2RotationOperator&                rotationOperator_;

   StokesGradientOperator   divT;
   StokesDivergenceOperator div;

   P1PSPGInvDiagOperator pspg_inv_diag_;

   StabilizationOperator_T stabOp_;

   SchurOperator_T schurOperator_;

   P2P1TaylorHoodFunction< real_t > tmp_;
   P2P1TaylorHoodFunction< real_t > tmpdst_;

   P2P1TaylorHoodFunction< idx_t > tmpAssembly_;
};

} // namespace hyteg