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

#include "hyteg/solvers/Solver.hpp"
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
   P2P1UzawaDampingFactorEstimationOperator( const std::shared_ptr< PrimitiveStorage >&                       storage,
                                             const std::shared_ptr< Solver< P2P1TaylorHoodStokesOperator > >& velocitySmoother,
                                             uint_t                                                           minLevel,
                                             uint_t                                                           maxLevel,
                                             const uint_t numGSIterationsVelocity = 2 )
   : Operator( storage, minLevel, maxLevel )
   , A( storage, minLevel, maxLevel )
   , mass_inv_diag_( storage, minLevel, maxLevel )
   , pspg_inv_diag_( storage, minLevel, maxLevel )
   , hasGlobalCells_( storage->hasGlobalCells() )
   , velocitySmoother_( velocitySmoother )
   , numGSIterationsVelocity_( numGSIterationsVelocity )
   , tmp_rhs_( "tmp_rhs_", storage, minLevel, maxLevel )
   , tmp_solution_( "tmp_solution_", storage, minLevel, maxLevel )
   , tmp_schur_( "tmp_schur", storage, minLevel, maxLevel )
   {}

   void apply( const P1Function< real_t >& src, const P1Function< real_t >& dst, const uint_t level, const DoFType flag ) const
   {
      tmp_solution_.u.interpolate( 0, level, All );
      tmp_solution_.v.interpolate( 0, level, All );
      tmp_solution_.w.interpolate( 0, level, All );

      A.divT_x.apply( src, tmp_rhs_.u, level, flag, Replace );
      A.divT_y.apply( src, tmp_rhs_.v, level, flag, Replace );
      if ( hasGlobalCells_ )
      {
         A.divT_z.apply( src, tmp_rhs_.w, level, flag, Replace );
      }

      for ( uint_t i = 0; i < numGSIterationsVelocity_; i++ )
      {
         velocitySmoother_->solve( A, tmp_solution_, tmp_rhs_, level );
      }

      A.div_x.apply( tmp_solution_.u, tmp_schur_, level, flag, Replace );
      A.div_y.apply( tmp_solution_.v, tmp_schur_, level, flag, Add );
      if ( hasGlobalCells_ )
      {
         A.div_z.apply( tmp_solution_.w, tmp_schur_, level, flag, Add );
      }

      mass_inv_diag_.apply( tmp_schur_, dst, level, flag, Replace );
   }

   P2P1TaylorHoodStokesOperator A;

   P1LumpedInvMassOperator                                   mass_inv_diag_;
   P1PSPGInvDiagOperator                                     pspg_inv_diag_;
   bool                                                      hasGlobalCells_;
   std::shared_ptr< Solver< P2P1TaylorHoodStokesOperator > > velocitySmoother_;

   uint_t numGSIterationsVelocity_;

   P2P1TaylorHoodFunction< real_t > tmp_rhs_;
   P2P1TaylorHoodFunction< real_t > tmp_solution_;

   P1Function< real_t > tmp_schur_;
};

} // namespace hyteg
