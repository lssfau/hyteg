/*
 * Copyright (c) 2017-2019 Nils Kohl.
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

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesBlockPreconditioner.hpp"
#include "hyteg/mixedoperators/P1ToP2Operator.hpp"
#include "hyteg/mixedoperators/P2ToP1Operator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"

namespace hyteg {

/// \brief Operator to wrap the application of
///
///   E = M^-1 * ( B * A^-1 * B^T )
///
/// where A^-1 corresponds to numGSIterationsVelocity symmetric or forward Gauss-Seidel iterations
/// and M^-1 is the inverse of the lumped mass-matrix.
///
/// A and B are the respective blocks of the P2-P1 discretized Stokes operator.
///
/// This wrapper is mainly used to approximate the relaxation parameter for the Uzawa smoother
/// by performing a few power iterations (solving the EV problem E v == lambda v).
///
/// See Drzisga et al: ON THE ANALYSIS OF BLOCK SMOOTHERS FOR SADDLE POINT PROBLEMS (sect. 5.2)
///
class P2P1UzawaDampingFactorEstimationOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
 public:
   P2P1UzawaDampingFactorEstimationOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                             uint_t                                     minLevel,
                                             uint_t                                     maxLevel,
                                             const bool                                 symmetricGSVelocity     = false,
                                             const uint_t                               numGSIterationsVelocity = 2 )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   , div_x( storage, minLevel, maxLevel )
   , div_y( storage, minLevel, maxLevel )
   , div_z( storage, minLevel, maxLevel )
   , divT_x( storage, minLevel, maxLevel )
   , divT_y( storage, minLevel, maxLevel )
   , divT_z( storage, minLevel, maxLevel )
   , mass_inv_diag_( storage, minLevel, maxLevel )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , symmetricGSVelocity_( symmetricGSVelocity )
   , numGSIterationsVelocity_( numGSIterationsVelocity )
   , tmp_rhs_u_( "tmp_rhs_u", storage, minLevel, maxLevel )
   , tmp_rhs_v_( "tmp_rhs_v", storage, minLevel, maxLevel )
   , tmp_rhs_w_( "tmp_rhs_w", storage, minLevel, maxLevel )
   , tmp_solution_u_( "tmp_solution_u", storage, minLevel, maxLevel )
   , tmp_solution_v_( "tmp_solution_v", storage, minLevel, maxLevel )
   , tmp_solution_w_( "tmp_solution_w", storage, minLevel, maxLevel )
   , tmp_schur_( "tmp_schur", storage, minLevel, maxLevel )
   {}

   void apply( const P1Function< real_t >& src,
               const P1Function< real_t >& dst,
               const uint_t                            level,
               const DoFType                           flag ) const
   {
      tmp_solution_u_.interpolate( 0, level, All );
      tmp_solution_v_.interpolate( 0, level, All );
      tmp_solution_w_.interpolate( 0, level, All );

      divT_x.apply( src, tmp_rhs_u_, level, flag, Replace );
      divT_y.apply( src, tmp_rhs_v_, level, flag, Replace );
      if ( hasGlobalCells_ )
      {
         divT_z.apply( src, tmp_rhs_w_, level, flag, Replace );
      }

      for ( uint_t i = 0; i < numGSIterationsVelocity_; i++ )
      {
         A.smooth_gs( tmp_solution_u_, tmp_rhs_u_, level, Inner );
         A.smooth_gs( tmp_solution_v_, tmp_rhs_v_, level, Inner );
         if ( hasGlobalCells_ )
         {
            A.smooth_gs( tmp_solution_w_, tmp_rhs_w_, level, Inner );
         }

         if ( symmetricGSVelocity_ )
         {
            A.smooth_gs_backwards( tmp_solution_u_, tmp_rhs_u_, level, Inner );
            A.smooth_gs_backwards( tmp_solution_v_, tmp_rhs_v_, level, Inner );
            if ( hasGlobalCells_ )
            {
               A.smooth_gs_backwards( tmp_solution_w_, tmp_rhs_w_, level, Inner );
            }
         }
      }

      div_x.apply( tmp_solution_u_, tmp_schur_, level, flag, Replace );
      div_y.apply( tmp_solution_v_, tmp_schur_, level, flag, Add );
      if ( hasGlobalCells_ )
      {
         div_z.apply( tmp_solution_w_, tmp_schur_, level, flag, Add );
      }

      mass_inv_diag_.apply( tmp_schur_, dst, level, flag, Replace );

   }

   P2ConstantLaplaceOperator   A;
   P2ToP1ConstantDivxOperator  div_x;
   P2ToP1ConstantDivyOperator  div_y;
   P2ToP1ConstantDivzOperator  div_z;
   P1ToP2ConstantDivTxOperator divT_x;
   P1ToP2ConstantDivTyOperator divT_y;
   P1ToP2ConstantDivTzOperator divT_z;

   /// this operator is need in the uzawa smoother
   P1PSPGInvDiagOperator pspg_inv_diag_;
   P1LumpedInvMassOperator mass_inv_diag_;
   bool                  hasGlobalCells_;
   bool symmetricGSVelocity_;
   uint_t numGSIterationsVelocity_;

   P2Function< real_t > tmp_rhs_u_;
   P2Function< real_t > tmp_rhs_v_;
   P2Function< real_t > tmp_rhs_w_;
   P2Function< real_t > tmp_solution_u_;
   P2Function< real_t > tmp_solution_v_;
   P2Function< real_t > tmp_solution_w_;
   P1Function< real_t > tmp_schur_;
};

} // namespace hyteg
