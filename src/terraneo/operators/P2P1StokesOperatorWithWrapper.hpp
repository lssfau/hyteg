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

#include "hyteg/p1functionspace/P1Petsc.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/p2functionspace/P2RotationOperator.hpp"
#include "hyteg_operators/operators/divergence/P2VectorToP1ElementwiseDivergenceCompressibleIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/divergence/P2VectorToP1ElementwiseDivergenceRotationIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/gradient/P1ToP2VectorElementwiseGradientRotationIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

namespace hyteg {

/***************************************************************************************************
NOTE: Stokes operator is wrapped with a FreeSlip Wrapper
      Changes the linear system from $[K] {u}           = {f} $
                                  to $[W] [K] [W]^T {u} = [W] {f}$
    
      where [W] can be a projection or a rotation operator
***************************************************************************************************/

template < typename FreeslipWrapperType,
           typename ViscosityFunctionType,
           typename DensityFunctionType,
           typename ViscosityInvFunctionType,
           typename StokesBaseOperatorType,
           typename SchurOperatorType,
           typename DivCompressibleOperatorType >
class P2P1FullStokesFreeslipWrapperTemplate
: public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   using FreeslipWrapper_T    = FreeslipWrapperType;
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

   class P2ViscousOperatorWithFreeslipWrapper : public Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >,
                                                public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >
   {
    public:
      P2ViscousOperatorWithFreeslipWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                                            uint_t                                     minLevel,
                                            uint_t                                     maxLevel,
                                            P2P1FullStokesFreeslipWrapperTemplate&     outerWrapper )
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
         outerWrapper_.freeslipWrapperOperator_.manipulate( outerWrapper_.tmp_->uvw(), level, FreeslipBoundary, true );

         {
            outerWrapper_.stokesOperator_.getA().apply(
                outerWrapper_.tmp_->uvw(), outerWrapper_.tmpdst_->uvw(), level, flag, updateType );
         }

         outerWrapper_.freeslipWrapperOperator_.manipulate( outerWrapper_.tmpdst_->uvw(), level, FreeslipBoundary );
         dst.assign( { real_c( 1.0 ) }, { outerWrapper_.tmpdst_->uvw() }, level, flag );
      }

      void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                     const P2VectorFunction< idx_t >&            src,
                     const P2VectorFunction< idx_t >&            dst,
                     size_t                                      level,
                     DoFType                                     flag ) const override
      {
         auto matProxyOp = mat->createCopy();
         outerWrapper_.stokesOperator_.getA().toMatrix( matProxyOp, src, dst, level, flag );

         const P2VectorFunction< idx_t >& wrapperSubstFunc = [&] {
            if constexpr ( std::is_same_v< FreeslipWrapper_T, P2ProjectNormalOperator > )
            {
               return src;
            }
            else if constexpr ( std::is_same_v< FreeslipWrapper_T, P2RotationOperator > )
            {
               return outerWrapper_.tmpAssembly_->uvw();
            }
            else
            {
               WALBERLA_ABORT( "Unknown type" );
            }
         }();

         auto matProxyProjectionPost = mat->createCopy();

         outerWrapper_.freeslipWrapperOperator_.toMatrix(
             matProxyProjectionPost, wrapperSubstFunc, level, FreeslipBoundary, false );

         std::vector< std::shared_ptr< SparseMatrixProxy > > matrices;
         matrices.push_back( matProxyProjectionPost );
         matrices.push_back( matProxyOp );

         //  if constexpr ( std::is_same_v< FreeslipWrapper_T, P2ProjectNormalOperator > )
         //  {
         auto matProxyProjectionPre = mat->createCopy();
         outerWrapper_.freeslipWrapperOperator_.toMatrix(
             matProxyProjectionPre, wrapperSubstFunc, level, FreeslipBoundary, true );

         matrices.push_back( matProxyProjectionPre );
         //  }

         mat->createFromMatrixProduct( matrices );
      }

      void computeInverseDiagonalOperatorValues() override { outerWrapper_.stokesOperator_.getA().computeInverseDiagonalOperatorValues(); }

      std::shared_ptr< P2VectorFunction< real_t > > getInverseDiagonalValues() const override
      {
         return outerWrapper_.stokesOperator_.getA().getInverseDiagonalValues();
      }

      P2P1FullStokesFreeslipWrapperTemplate& outerWrapper_;
   };

   class P1ToP2GradientOperatorWithFreeslipWrapper : public Operator< P1Function< real_t >, P2VectorFunction< real_t > >
   {
    public:
      P1ToP2GradientOperatorWithFreeslipWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                                                 uint_t                                     minLevel,
                                                 uint_t                                     maxLevel,
                                                 P2P1FullStokesFreeslipWrapperTemplate&     outerWrapper )
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

         outerWrapper_.freeslipWrapperOperator_.manipulate( outerWrapper_.tmp_->uvw(), level, FreeslipBoundary );
         dst.assign( { updateType == Add ? real_c(1.0) : real_c(0.0), real_c(1.0) }, { dst, outerWrapper_.tmp_->uvw() }, level, flag );
      }

      P2P1FullStokesFreeslipWrapperTemplate& outerWrapper_;
   };

   class P2ToP1DivergenceOperatorWithFreeslipWrapper : public Operator< P2VectorFunction< real_t >, P1Function< real_t > >
   {
    public:
      P2ToP1DivergenceOperatorWithFreeslipWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                                                   uint_t                                     minLevel,
                                                   uint_t                                     maxLevel,
                                                   P2P1FullStokesFreeslipWrapperTemplate&     outerWrapper )
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
         outerWrapper_.freeslipWrapperOperator_.manipulate( outerWrapper_.tmp_->uvw(), level, FreeslipBoundary, true );

         if ( !( outerWrapper_.frozenVelocity_ ) )
         {
            outerWrapper_.divCompressibleOperator_->apply( outerWrapper_.tmp_->uvw(), dst, level, flag, updateType );
         }
         else
         {
            outerWrapper_.stokesOperator_.getB().apply( outerWrapper_.tmp_->uvw(), dst, level, flag, updateType );
         }
      }

      P2P1FullStokesFreeslipWrapperTemplate& outerWrapper_;
   };

   using ViscousOperatorWrapped_T = P2ViscousOperatorWithFreeslipWrapper;
   using GradOperatorWrapped_T    = P1ToP2GradientOperatorWithFreeslipWrapper;
   using DivOperatorWrapped_T     = P2ToP1DivergenceOperatorWithFreeslipWrapper;

   using VelocityOperator_T  = ViscousOperatorWrapped_T;
   using ViscousOperatorFS_T = ViscousOperatorWrapped_T;

   P2P1FullStokesFreeslipWrapperTemplate( const std::shared_ptr< PrimitiveStorage >&                storage,
                                          uint_t                                                    minLevel,
                                          uint_t                                                    maxLevel,
                                          const ViscosityFunction_T&                                mu,
                                          const ViscosityInvFunction_T&                             muInv,
                                          FreeslipWrapper_T&                                        freeslipWrapperOperator,
                                          BoundaryCondition                                         bcVelocity,
                                          const DensityFunction_T&                                  rhoFE,
                                          const bool                                                frozenVelocity = true,
                                          const std::shared_ptr< P2P1TaylorHoodFunction< real_t > > tmp            = nullptr,
                                          const std::shared_ptr< P2P1TaylorHoodFunction< real_t > > tmpdst         = nullptr,
                                          const std::shared_ptr< P2P1TaylorHoodFunction< idx_t > >  tmpAssembly    = nullptr )
   : Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >( storage, minLevel, maxLevel )
   , frozenVelocity_( frozenVelocity )
   , stokesOperator_( storage, minLevel, maxLevel, mu )
   , schurOperator_( storage, minLevel, maxLevel, muInv )
   , freeslipWrapperOperator_( freeslipWrapperOperator )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , massOperator( storage, minLevel, maxLevel )
   , gradRotationOperator_( storage, minLevel, maxLevel, *this )
   , divRotationOperator_( storage, minLevel, maxLevel, *this )
   , viscousOperatorWrapped_( storage, minLevel, maxLevel, *this )
   {
      if ( tmp == nullptr )
      {
         tmp_ = std::make_shared< P2P1TaylorHoodFunction< real_t > >(
             "tmp__P2P1StokesFullIcosahedralShellMapOperatorFS", storage, minLevel, maxLevel, bcVelocity );
      }
      else
      {
         tmp_ = tmp;
      }

      if ( tmpdst == nullptr )
      {
         tmpdst_ = std::make_shared< P2P1TaylorHoodFunction< real_t > >(
             "tmpdst__P2P1StokesFullIcosahedralShellMapOperatorFS", storage, minLevel, maxLevel, bcVelocity );
      }
      else
      {
         tmpdst_ = tmpdst;
      }

      if ( tmpAssembly == nullptr )
      {
         tmpAssembly_ = std::make_shared< P2P1TaylorHoodFunction< idx_t > >(
             "tmpAssembly__P2P1StokesFullIcosahedralShellMapOperatorFS", storage, minLevel, maxLevel, bcVelocity );
      }
      else
      {
         tmpAssembly_ = tmpAssembly;
      }

      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         tmpAssembly_->enumerate( level );
      }

      if ( !frozenVelocity_ )
      {
         divCompressibleOperator_ = std::make_shared< DivCompressibleOperator_T >( storage, minLevel, maxLevel, rhoFE );
      }

      stokesOperator_.getA().computeInverseDiagonalOperatorValues();
      schurOperator_.computeInverseDiagonalOperatorValues();
   }

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag,
               const UpdateType                        updateType = Replace ) const override
   {
      tmp_->assign( { real_c( 1.0 ) }, { src }, level, All );
      freeslipWrapperOperator_.manipulate( tmp_->uvw(), level, FreeslipBoundary, true );

      {
         stokesOperator_.getA().apply( tmp_->uvw(), tmpdst_->uvw(), level, flag, updateType );
         if ( !frozenVelocity_ )
         {
            divCompressibleOperator_->apply( tmp_->uvw(), tmpdst_->p(), level, flag, updateType );
         }
         else
         {
            stokesOperator_.getB().apply( tmp_->uvw(), tmpdst_->p(), level, flag, updateType );
         }
         stokesOperator_.getBT().apply( tmp_->p(), tmpdst_->uvw(), level, flag, Add );
      }

      freeslipWrapperOperator_.manipulate( tmpdst_->uvw(), level, FreeslipBoundary );
      dst.assign( { real_c( 1.0 ) }, { *tmpdst_ }, level, flag );

      vertexdof::projectMean( dst.p(), level );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2P1TaylorHoodFunction< idx_t >&      src,
                  const P2P1TaylorHoodFunction< idx_t >&      dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      auto matProxyOp = mat->createCopy();

      {
         stokesOperator_.getA().toMatrix( matProxyOp, src.uvw(), dst.uvw(), level, flag );
         if ( !frozenVelocity_ )
         {
            divCompressibleOperator_->toMatrix( matProxyOp, src.uvw(), dst.p(), level, flag );
         }
         else
         {
            stokesOperator_.getB().toMatrix( matProxyOp, src.uvw(), dst.p(), level, flag );
         }
         stokesOperator_.getBT().toMatrix( matProxyOp, src.p(), dst.uvw(), level, flag );
      }

      auto matProxyProjectionPost = mat->createCopy();
      freeslipWrapperOperator_.toMatrix( matProxyProjectionPost, tmpAssembly_->uvw(), level, FreeslipBoundary, false );
      saveIdentityOperator( tmpAssembly_->p(), matProxyProjectionPost, level, All );

      std::vector< std::shared_ptr< SparseMatrixProxy > > matrices;
      matrices.push_back( matProxyProjectionPost );
      matrices.push_back( matProxyOp );

      auto matProxyProjectionPre = mat->createCopy();
      freeslipWrapperOperator_.toMatrix( matProxyProjectionPre, tmpAssembly_->uvw(), level, FreeslipBoundary, true );
      saveIdentityOperator( tmpAssembly_->p(), matProxyProjectionPre, level, All );
      matrices.push_back( matProxyProjectionPre );

      mat->createFromMatrixProduct( matrices );
   }

   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > tmp_;
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > tmpdst_;
   std::shared_ptr< P2P1TaylorHoodFunction< idx_t > >  tmpAssembly_;

   const bool frozenVelocity_;

   StokesOperator_T   stokesOperator_;
   SchurOperator_T    schurOperator_;
   FreeslipWrapper_T& freeslipWrapperOperator_;

   std::shared_ptr< DivCompressibleOperator_T > divCompressibleOperator_;

   P1PSPGInvDiagOperator pspg_inv_diag_;

   P2ElementwiseBlendingMassOperator massOperator;

   GradOperatorWrapped_T    gradRotationOperator_;
   DivOperatorWrapped_T     divRotationOperator_;
   ViscousOperatorWrapped_T viscousOperatorWrapped_;

   const ViscousOperatorWrapped_T& getA() const { return viscousOperatorWrapped_; }
   const DivOperatorWrapped_T&     getB() const { return divRotationOperator_; }
   const GradOperatorWrapped_T&    getBT() const { return gradRotationOperator_; }

   const SchurOperator_T& getSchur() const { return schurOperator_; }
   const StabOperator_T&  getStab() const { return stokesOperator_.getStab(); }

   const GradOperator_T divT = stokesOperator_.getBT();
   const DivOperator_T  div  = stokesOperator_.getB();

   ViscousOperatorWrapped_T& getA() { return viscousOperatorWrapped_; }
};

using P2P1StokesFullIcosahedralShellMapOperatorWithProjection =
    P2P1FullStokesFreeslipWrapperTemplate< P2ProjectNormalOperator,
                                           P2Function< real_t >,
                                           P2Function< real_t >,
                                           P1Function< real_t >,
                                           operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator,
                                           operatorgeneration::P1ElementwiseKMassIcosahedralShellMap,
                                           operatorgeneration::P2VectorToP1ElementwiseDivergenceCompressibleIcosahedralShellMap >;

using P2P1StokesP1ViscosityFullIcosahedralShellMapOperatorWithProjection =
    P2P1FullStokesFreeslipWrapperTemplate< P2ProjectNormalOperator,
                                           P1Function< real_t >,
                                           P2Function< real_t >,
                                           P1Function< real_t >,
                                           operatorgeneration::P2P1StokesFullP1ViscosityIcosahedralShellMapOperator,
                                           operatorgeneration::P1ElementwiseKMassIcosahedralShellMap,
                                           operatorgeneration::P2VectorToP1ElementwiseDivergenceCompressibleIcosahedralShellMap >;

using P2P1StokesP0ViscosityFullIcosahedralShellMapOperatorWithProjection =
    P2P1FullStokesFreeslipWrapperTemplate< P2ProjectNormalOperator,
                                           P0Function< real_t >,
                                           P2Function< real_t >,
                                           P1Function< real_t >,
                                           operatorgeneration::P2P1StokesFullP0ViscosityIcosahedralShellMapOperator,
                                           operatorgeneration::P1ElementwiseKMassIcosahedralShellMap,
                                           operatorgeneration::P2VectorToP1ElementwiseDivergenceCompressibleIcosahedralShellMap >;

using P2P1StokesFullIcosahedralShellMapOperatorWithRotation =
    P2P1FullStokesFreeslipWrapperTemplate< P2RotationOperator,
                                           P2Function< real_t >,
                                           P2Function< real_t >,
                                           P1Function< real_t >,
                                           operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator,
                                           operatorgeneration::P1ElementwiseKMassIcosahedralShellMap,
                                           operatorgeneration::P2VectorToP1ElementwiseDivergenceCompressibleIcosahedralShellMap >;

using P2P1StokesP1ViscosityFullIcosahedralShellMapOperatorWithRotation =
    P2P1FullStokesFreeslipWrapperTemplate< P2RotationOperator,
                                           P1Function< real_t >,
                                           P2Function< real_t >,
                                           P1Function< real_t >,
                                           operatorgeneration::P2P1StokesFullP1ViscosityIcosahedralShellMapOperator,
                                           operatorgeneration::P1ElementwiseKMassIcosahedralShellMap,
                                           operatorgeneration::P2VectorToP1ElementwiseDivergenceCompressibleIcosahedralShellMap >;

using P2P1StokesP0ViscosityFullIcosahedralShellMapOperatorWithRotation =
    P2P1FullStokesFreeslipWrapperTemplate< P2RotationOperator,
                                           P0Function< real_t >,
                                           P2Function< real_t >,
                                           P1Function< real_t >,
                                           operatorgeneration::P2P1StokesFullP0ViscosityIcosahedralShellMapOperator,
                                           operatorgeneration::P1ElementwiseKMassIcosahedralShellMap,
                                           operatorgeneration::P2VectorToP1ElementwiseDivergenceCompressibleIcosahedralShellMap >;
} // namespace hyteg