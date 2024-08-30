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

#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesEpsilonOperator.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

#include "terraneo/operators/P2StokesABlockWithProjection.hpp"

namespace hyteg {

/***************************************************************************************************
NOTE: Here FS denotes FreeSlip, Normal Stokes operator is wrapped with a FreeSlip Projection Wrapper
      Changes the linear system from $Ku    = f $
                                  to $PKP^Tu = Pf$
***************************************************************************************************/
class P2P1StokesFullIcosahedralShellMapOperatorFS
: public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator::ViscousOperator_T       ViscousOperator_T;
   typedef P2ViscousIcosahedralShellMapOperatorFS                                                 ViscousOperatorFS_T;
   typedef ViscousOperatorFS_T                                                                    VelocityOperator_T;
   typedef operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator::DivergenceOperator_T    DivOperator_T;
   typedef operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator::GradientOperator_T      GradOperator_T;
   typedef operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator::StabilizationOperator_T StabOperator_T;

   typedef operatorgeneration::P1ElementwiseKMassIcosahedralShellMap SchurOperator_T;

   P2P1StokesFullIcosahedralShellMapOperatorFS( const std::shared_ptr< PrimitiveStorage >& storage,
                                                uint_t                                     minLevel,
                                                uint_t                                     maxLevel,
                                                const P2Function< real_t >&                mu,
                                                const P1Function< real_t >&                muInv,
                                                P2ProjectNormalOperator&                   projectNormal,
                                                BoundaryCondition                          bcVelocity )
   : Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >( storage, minLevel, maxLevel )
   , tmp_( "tmp__P2P1StokesFullIcosahedralShellMapOperatorFS", storage, minLevel, maxLevel, bcVelocity )
   , tmp_Vec( "tmp_Vec_P2P1StokesFullIcosahedralShellMapOperatorFS", storage, minLevel, maxLevel )
   , tmp_P2( "tmp_P2_P2P1StokesFullIcosahedralShellMapOperatorFS", storage, minLevel, maxLevel )
   , StokesOp( storage, minLevel, maxLevel, mu )
   , schurOperator( storage, minLevel, maxLevel, muInv )
   , projectNormal_( projectNormal )
   , viscousFSOp( storage, minLevel, maxLevel, mu, projectNormal_, bcVelocity )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , massOperator( storage, minLevel, maxLevel )
   {
      StokesOp.getA().computeInverseDiagonalOperatorValues();
      viscousFSOp.computeInverseDiagonalOperatorValues();
      schurOperator.computeInverseDiagonalOperatorValues();
   }

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag,
               const UpdateType                        updateType = Replace ) const
   {
      // hyteg::removeRotationalModes( massOperator, src.uvw(), tmp_Vec, tmp_P2, level );
      // vertexdof::projectMean( src.p(), level );

      tmp_.assign( { 1 }, { src }, level, All );
      projectNormal_.project( tmp_, level, FreeslipBoundary );

      StokesOp.apply( tmp_, dst, level, flag, updateType );
      projectNormal_.project( dst, level, FreeslipBoundary );

      vertexdof::projectMean( dst.p(), level );
   }

   P2P1TaylorHoodFunction< real_t > tmp_;

   P2VectorFunction< real_t > tmp_Vec;
   P2Function< real_t >       tmp_P2;

   operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator StokesOp;
   SchurOperator_T                                               schurOperator;
   P2ProjectNormalOperator&                                      projectNormal_;
   ViscousOperatorFS_T                                           viscousFSOp;

   P1PSPGInvDiagOperator pspg_inv_diag_;

   P2ElementwiseBlendingMassOperator massOperator;

   const ViscousOperatorFS_T& getA() const { return viscousFSOp; }
   const DivOperator_T&       getB() const { return StokesOp.getB(); }
   const GradOperator_T&      getBT() const { return StokesOp.getBT(); }
   const SchurOperator_T&     getSchur() const { return schurOperator; }
   const StabOperator_T&      getStab() const { return StokesOp.getStab(); }

   const GradOperator_T divT = StokesOp.getBT();
   const DivOperator_T  div  = StokesOp.getB();

   ViscousOperatorFS_T& getA() { return viscousFSOp; }
};

class P2P1StokesP1ViscosityFullIcosahedralShellMapOperatorFS
: public Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef operatorgeneration::P2P1StokesFullP1ViscosityIcosahedralShellMapOperator    StokesOperatorBase_T;
   // typedef operatorgeneration::P2P1StokesEpsilonP1ViscosityIcosahedralShellMapOperator StokesOperatorBase_T;
   typedef StokesOperatorBase_T                                                        StokesOperator_T;
   typedef StokesOperatorBase_T::ViscousOperator_T                                     ViscousOperator_T;
   typedef P2ABlockP1ViscousIcosahedralShellMapOperatorFS                              ViscousOperatorFS_T;
   typedef ViscousOperatorFS_T                                                         VelocityOperator_T;
   typedef StokesOperatorBase_T::DivergenceOperator_T                                  DivOperator_T;
   typedef StokesOperatorBase_T::GradientOperator_T                                    GradOperator_T;
   typedef StokesOperatorBase_T::StabilizationOperator_T                               StabOperator_T;

   typedef operatorgeneration::P1ElementwiseKMassIcosahedralShellMap SchurOperator_T;

   P2P1StokesP1ViscosityFullIcosahedralShellMapOperatorFS( const std::shared_ptr< PrimitiveStorage >& storage,
                                                           uint_t                                     minLevel,
                                                           uint_t                                     maxLevel,
                                                           const P1Function< real_t >&                mu,
                                                           const P1Function< real_t >&                muInv,
                                                           P2ProjectNormalOperator&                   projectNormal,
                                                           BoundaryCondition                          bcVelocity )
   : Operator< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodFunction< real_t > >( storage, minLevel, maxLevel )
   , tmp_( "tmp__P2P1StokesFullIcosahedralShellMapOperatorFS", storage, minLevel, maxLevel, bcVelocity )
   , tmp_Vec( "tmp_Vec_P2P1StokesFullIcosahedralShellMapOperatorFS", storage, minLevel, maxLevel )
   , tmp_P2( "tmp_P2_P2P1StokesFullIcosahedralShellMapOperatorFS", storage, minLevel, maxLevel )
   , StokesOp( storage, minLevel, maxLevel, mu )
   , schurOperator( storage, minLevel, maxLevel, muInv )
   , projectNormal_( projectNormal )
   , viscousFSOp( storage, minLevel, maxLevel, mu, projectNormal_, bcVelocity )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , massOperator( storage, minLevel, maxLevel )
   {
      StokesOp.getA().computeInverseDiagonalOperatorValues();
      viscousFSOp.computeInverseDiagonalOperatorValues();
      schurOperator.computeInverseDiagonalOperatorValues();
   }

   void apply( const P2P1TaylorHoodFunction< real_t >& src,
               const P2P1TaylorHoodFunction< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag,
               const UpdateType                        updateType = Replace ) const
   {
      // hyteg::removeRotationalModes( massOperator, src.uvw(), tmp_Vec, tmp_P2, level );
      // vertexdof::projectMean( src.p(), level );

      tmp_.assign( { 1 }, { src }, level, All );
      projectNormal_.project( tmp_, level, FreeslipBoundary );

      StokesOp.apply( tmp_, dst, level, flag, updateType );
      projectNormal_.project( dst, level, FreeslipBoundary );

      vertexdof::projectMean( dst.p(), level );
   }

   P2P1TaylorHoodFunction< real_t > tmp_;

   P2VectorFunction< real_t > tmp_Vec;
   P2Function< real_t >       tmp_P2;

   StokesOperator_T         StokesOp;
   SchurOperator_T          schurOperator;
   P2ProjectNormalOperator& projectNormal_;
   ViscousOperatorFS_T      viscousFSOp;

   P1PSPGInvDiagOperator pspg_inv_diag_;

   P2ElementwiseBlendingMassOperator massOperator;

   const ViscousOperatorFS_T& getA() const { return viscousFSOp; }
   const DivOperator_T&       getB() const { return StokesOp.getB(); }
   const GradOperator_T&      getBT() const { return StokesOp.getBT(); }
   const SchurOperator_T&     getSchur() const { return schurOperator; }
   const StabOperator_T&      getStab() const { return StokesOp.getStab(); }

   const GradOperator_T divT = StokesOp.getBT();
   const DivOperator_T  div  = StokesOp.getB();

   ViscousOperatorFS_T& getA() { return viscousFSOp; }
};
} // namespace hyteg