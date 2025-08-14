/*
 * Copyright (c) 2024 Ponsuganth Ilangovan P, Andreas Burkhart.
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

#include "hyteg/p1functionspace/P1Petsc.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg_operators/operators/divergence/P2VectorToP1ElementwiseDivergenceCompressibleIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesEpsilonOperator.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

// PETSc
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"

namespace hyteg {

/***************************************************************************************************
NOTE: Here FS denotes FreeSlip, Normal Stokes operator is wrapped with a FreeSlip Projection Wrapper
      Changes the linear system from $Ku    = f $
                                  to $PKP^Tu = Pf$
***************************************************************************************************/

template < typename ViscosityFunctionType,
           typename DensityFunctionType,
           typename ViscosityInvFunctionType,
           typename StokesBaseOperatorType,
           typename SchurOperatorType,
           typename DivCompressibleOperatorType >
class P2P1FullStokesProjectionFSTemplate : public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   using StokesOperatorBase_T = StokesBaseOperatorType;
   using SchurOperator_T      = SchurOperatorType;

   using ViscosityFunction_T    = ViscosityFunctionType;
   using DensityFunction_T      = DensityFunctionType;
   using ViscosityInvFunction_T = ViscosityInvFunctionType;

   using StokesOperator_T = StokesOperatorBase_T;

   using ViscousOperator_T = typename StokesOperatorBase_T::ViscousOperator_T;
   using DivOperator_T     = typename StokesOperatorBase_T::DivergenceOperator_T;
   using GradOperator_T    = typename StokesOperatorBase_T::GradientOperator_T;
   using StabOperator_T    = typename StokesOperatorBase_T::StabilizationOperator_T;

   using DivCompressibleOperator_T = DivCompressibleOperatorType;

   class P2ViscousOperatorWithProjection : public Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >,
                                           public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >
   {
    public:
      P2ViscousOperatorWithProjection( const std::shared_ptr< PrimitiveStorage >& storage,
                                       uint_t                                     minLevel,
                                       uint_t                                     maxLevel,
                                       P2P1FullStokesProjectionFSTemplate&        outerWrapper )
      : Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >( storage, minLevel, maxLevel )
      , outerWrapper_( outerWrapper )
      {}

      void apply( const P2VectorFunction< real_t >& src,
                  const P2VectorFunction< real_t >& dst,
                  const uint_t                      level,
                  const DoFType                     flag,
                  const UpdateType                  updateType = Replace ) const override
      {
         outerWrapper_.tmp_->uvw().assign( { real_c( 1.0 ) }, { src }, level, All );

         outerWrapper_.projectNormal_.project( outerWrapper_.tmp_->uvw(), level, FreeslipBoundary );
         outerWrapper_.stokesOperator_.getA().apply( outerWrapper_.tmp_->uvw(), dst, level, flag, updateType );
         outerWrapper_.projectNormal_.project( dst, level, FreeslipBoundary );
      }

      void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                     const P2VectorFunction< idx_t >&            src,
                     const P2VectorFunction< idx_t >&            dst,
                     size_t                                      level,
                     DoFType                                     flag ) const override
      {
         auto matProxyOp = mat->createCopy();
         outerWrapper_.stokesOperator_.getA().toMatrix( matProxyOp, src, dst, level, flag );

         auto matProxyProjectionPre  = mat->createCopy();
         auto matProxyProjectionPost = mat->createCopy();

         outerWrapper_.projectNormal_.toMatrix( matProxyProjectionPost, src, level, FreeslipBoundary );
         outerWrapper_.projectNormal_.toMatrix( matProxyProjectionPre, src, level, FreeslipBoundary );

         std::vector< std::shared_ptr< SparseMatrixProxy > > matrices;

         matrices.push_back( matProxyProjectionPost );
         matrices.push_back( matProxyOp );
         matrices.push_back( matProxyProjectionPre );

         mat->createFromMatrixProduct( matrices );
      }

      void computeInverseDiagonalOperatorValues() { outerWrapper_.stokesOperator_.getA().computeInverseDiagonalOperatorValues(); }

      std::shared_ptr< P2VectorFunction< real_t > > getInverseDiagonalValues() const
      {
         return outerWrapper_.stokesOperator_.getA().getInverseDiagonalValues();
      }

      P2P1FullStokesProjectionFSTemplate& outerWrapper_;
   };

   class P1ToP2GradientOperatorWithProjection : public Operator< P1Function< real_t >, P2VectorFunction< real_t > >
   {
    public:
      P1ToP2GradientOperatorWithProjection( const std::shared_ptr< PrimitiveStorage >& storage,
                                            uint_t                                     minLevel,
                                            uint_t                                     maxLevel,
                                            P2P1FullStokesProjectionFSTemplate&        outerWrapper )
      : Operator< P1Function< real_t >, P2VectorFunction< real_t > >( storage, minLevel, maxLevel )
      , outerWrapper_( outerWrapper )
      {}

      void apply( const P1Function< real_t >&       src,
                  const P2VectorFunction< real_t >& dst,
                  const uint_t                      level,
                  const DoFType                     flag,
                  const UpdateType                  updateType = Replace ) const override
      {
         {
            outerWrapper_.stokesOperator_.getBT().apply( src, outerWrapper_.tmp_->uvw(), level, flag );
         }

         outerWrapper_.projectNormal_.project( outerWrapper_.tmp_->uvw(), level, FreeslipBoundary );
         dst.assign( { updateType == Add ? 1.0 : 0.0, 1.0 }, { dst, outerWrapper_.tmp_->uvw() }, level, flag );
      }

      P2P1FullStokesProjectionFSTemplate& outerWrapper_;
   };

   class P2ToP1DivergenceOperatorWithProjection : public Operator< P2VectorFunction< real_t >, P1Function< real_t > >
   {
    public:
      P2ToP1DivergenceOperatorWithProjection( const std::shared_ptr< PrimitiveStorage >& storage,
                                              uint_t                                     minLevel,
                                              uint_t                                     maxLevel,
                                              P2P1FullStokesProjectionFSTemplate&        outerWrapper )
      : Operator< P2VectorFunction< real_t >, P1Function< real_t > >( storage, minLevel, maxLevel )
      , outerWrapper_( outerWrapper )
      {}

      void apply( const P2VectorFunction< real_t >& src,
                  const P1Function< real_t >&       dst,
                  const uint_t                      level,
                  const DoFType                     flag,
                  const UpdateType                  updateType = Replace ) const override
      {
         outerWrapper_.tmp_->uvw().assign( { real_c( 1.0 ) }, { src }, level, All );
         outerWrapper_.projectNormal_.project( outerWrapper_.tmp_->uvw(), level, FreeslipBoundary );

         if ( !outerWrapper_.frozenVelocity_ )
         {
            outerWrapper_.divCompressibleOperator_->apply( outerWrapper_.tmp_->uvw(), dst, level, flag, updateType );
         }
         else
         {
            outerWrapper_.stokesOperator_.getB().apply( outerWrapper_.tmp_->uvw(), dst, level, flag, updateType );
         }
      }

      P2P1FullStokesProjectionFSTemplate& outerWrapper_;
   };

   using ViscousOperatorFS_T = P2ViscousOperatorWithProjection;
   using VelocityOperator_T  = ViscousOperatorFS_T;

   using GradientProjectionOperator_T   = P1ToP2GradientOperatorWithProjection;
   using DivergenceProjectionOperator_T = P2ToP1DivergenceOperatorWithProjection;

   P2P1FullStokesProjectionFSTemplate( const std::shared_ptr< PrimitiveStorage >&          storage,
                                       uint_t                                              minLevel,
                                       uint_t                                              maxLevel,
                                       const ViscosityFunction_T&                          mu,
                                       const ViscosityInvFunction_T&                       muInv,
                                       P2ProjectNormalOperator&                            projectNormal,
                                       BoundaryCondition                                   bcVelocity,
                                       const DensityFunction_T&                            rhoFE,
                                       const bool                                          frozenVelocity = true,
                                       std::shared_ptr< P2P1TaylorHoodFunction< real_t > > tmp            = nullptr,
                                       const real_t                                        rotFactor      = 0.0 )
   : Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >( storage, minLevel, maxLevel )
   , frozenVelocity_( frozenVelocity )
   , stokesOperator_( storage, minLevel, maxLevel, mu )
   , schurOperator_( storage, minLevel, maxLevel, muInv )
   , projectNormal_( projectNormal )
   , viscousOperatorWrapped_( storage, minLevel, maxLevel, *this )
   , gradOperatorWrapped_( storage, minLevel, maxLevel, *this )
   , divOperatorWrapped_( storage, minLevel, maxLevel, *this )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , massOperator( storage, minLevel, maxLevel )
   {
      if ( tmp == nullptr )
      {
         tmp_ = std::make_shared< P2P1TaylorHoodFunction< real_t > >(
             "tmp__P2P1FullStokesProjectionFSTemplate", storage, minLevel, maxLevel, bcVelocity );
      }
      else
      {
         tmp_ = tmp;
      }

      if ( !frozenVelocity_ )
      {
         divCompressibleOperator_ = std::make_shared< DivCompressibleOperator_T >( storage, minLevel, maxLevel, rhoFE );
      }

      stokesOperator_.getA().computeInverseDiagonalOperatorValues();
      viscousOperatorWrapped_.computeInverseDiagonalOperatorValues();
      schurOperator_.computeInverseDiagonalOperatorValues();
   }

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag,
               const UpdateType                        updateType = Replace ) const override
   {
      tmp_->assign( { 1 }, { src }, level, All );
      projectNormal_.project( *tmp_, level, FreeslipBoundary );

      stokesOperator_.getA().apply( tmp_->uvw(), dst.uvw(), level, flag, updateType );
      if ( !frozenVelocity_ )
      {
         divCompressibleOperator_->apply( tmp_->uvw(), dst.p(), level, flag, updateType );
      }
      else
      {
         stokesOperator_.getB().apply( tmp_->uvw(), dst.p(), level, flag, updateType );
      }
      stokesOperator_.getBT().apply( tmp_->p(), dst.uvw(), level, flag, Add );

      projectNormal_.project( dst, level, FreeslipBoundary );
      vertexdof::projectMean( dst.p(), level );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< idx_t >&      src,
                  const P2P1TaylorHoodFunction< idx_t >&      dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      auto matProxyOp = mat->createCopy();
      stokesOperator_.getA().toMatrix( matProxyOp, src.uvw(), dst.uvw(), level, flag );
      divT.toMatrix( matProxyOp, src.p(), dst.uvw(), level, flag );

      if ( !frozenVelocity_ )
      {
         divCompressibleOperator_->toMatrix( matProxyOp, src.uvw(), dst.p(), level, flag );
      }
      else
      {
         stokesOperator_.getB().toMatrix( matProxyOp, src.uvw(), dst.p(), level, flag );
      }

      auto matProxyProjectionPre  = mat->createCopy();
      auto matProxyProjectionPost = mat->createCopy();

      projectNormal_.toMatrix( matProxyProjectionPre, src.uvw(), dst.uvw(), level, FreeslipBoundary );
      projectNormal_.toMatrix( matProxyProjectionPost, src.uvw(), dst.uvw(), level, FreeslipBoundary );

      // save ID in pressure block
      saveIdentityOperator( dst.p(), matProxyProjectionPost, level, All );

      std::vector< std::shared_ptr< SparseMatrixProxy > > matrices;

      matrices.push_back( matProxyProjectionPost );
      matrices.push_back( matProxyOp );
      matrices.push_back( matProxyProjectionPre );

      mat->createFromMatrixProduct( matrices );
   }

   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > tmp_;

   const bool frozenVelocity_;

   StokesOperator_T         stokesOperator_;
   SchurOperator_T          schurOperator_;
   P2ProjectNormalOperator& projectNormal_;

   ViscousOperatorFS_T            viscousOperatorWrapped_;
   GradientProjectionOperator_T   gradOperatorWrapped_;
   DivergenceProjectionOperator_T divOperatorWrapped_;

   std::shared_ptr< DivCompressibleOperator_T > divCompressibleOperator_;

   P1PSPGInvDiagOperator pspg_inv_diag_;

   P2ElementwiseBlendingMassOperator massOperator;

   const ViscousOperatorFS_T&            getA() const { return viscousOperatorWrapped_; }
   const DivergenceProjectionOperator_T& getB() const { return divOperatorWrapped_; }
   const GradientProjectionOperator_T&   getBT() const { return gradOperatorWrapped_; }
   const SchurOperator_T&                getSchur() const { return schurOperator_; }
   const StabOperator_T&                 getStab() const { return stokesOperator_.getStab(); }

   const GradOperator_T divT = stokesOperator_.getBT();
   const DivOperator_T  div  = stokesOperator_.getB();

   ViscousOperatorFS_T& getA() { return viscousOperatorWrapped_; }
};

using P2P1StokesFullIcosahedralShellMapOperatorFS =
    P2P1FullStokesProjectionFSTemplate< P2Function< real_t >,
                                        P2Function< real_t >,
                                        P1Function< real_t >,
                                        operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator,
                                        operatorgeneration::P1ElementwiseKMassIcosahedralShellMap,
                                        operatorgeneration::P2VectorToP1ElementwiseDivergenceCompressibleIcosahedralShellMap >;

using P2P1StokesP1ViscosityFullIcosahedralShellMapOperatorFS =
    P2P1FullStokesProjectionFSTemplate< P1Function< real_t >,
                                        P2Function< real_t >,
                                        P1Function< real_t >,
                                        operatorgeneration::P2P1StokesFullP1ViscosityIcosahedralShellMapOperator,
                                        operatorgeneration::P1ElementwiseKMassIcosahedralShellMap,
                                        operatorgeneration::P2VectorToP1ElementwiseDivergenceCompressibleIcosahedralShellMap >;

using P2P1StokesP0ViscosityFullIcosahedralShellMapOperatorFS =
    P2P1FullStokesProjectionFSTemplate< P0Function< real_t >,
                                        P2Function< real_t >,
                                        P1Function< real_t >,
                                        operatorgeneration::P2P1StokesFullP0ViscosityIcosahedralShellMapOperator,
                                        operatorgeneration::P1ElementwiseKMassIcosahedralShellMap,
                                        operatorgeneration::P2VectorToP1ElementwiseDivergenceCompressibleIcosahedralShellMap >;

} // namespace hyteg