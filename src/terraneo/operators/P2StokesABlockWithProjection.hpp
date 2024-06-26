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
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

namespace hyteg {

/***************************************************************************************************
NOTE: Here FS denotes FreeSlip, Stokes A block operator is wrapped with a FreeSlip Projection Wrapper
      Changes the linear system from $Ku    = f $
                                  to $PKP^T = Pf$
***************************************************************************************************/
class P2ViscousIcosahedralShellMapOperatorFS : public Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >
{
 public:
   typedef operatorgeneration::P2P1StokesFullIcosahedralShellMapOperator::ViscousOperator_T ViscousOperator_T;

   P2ViscousIcosahedralShellMapOperatorFS( const std::shared_ptr< PrimitiveStorage >& storage,
                                           uint_t                                     minLevel,
                                           uint_t                                     maxLevel,
                                           const P2Function< real_t >&                mu,
                                           P2ProjectNormalOperator&                   projectNormal,
                                           BoundaryCondition                          bcVelocity )
   : Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >( storage, minLevel, maxLevel )
   , tmp_( "tmp__P2ViscousIcosahedralShellMapOperatorFS", storage, minLevel, maxLevel, bcVelocity )
   , viscousOperator( storage, minLevel, maxLevel, mu )
   , projectNormal_( projectNormal )
   {}

   void apply( const P2VectorFunction< real_t >& src,
               const P2VectorFunction< real_t >& dst,
               const uint_t                      level,
               const DoFType                     flag,
               const UpdateType                  updateType = Replace ) const
   {
      tmp_.assign( { 1 }, { src }, level, All );
      // projectNormal_.project( tmp_, level, FreeslipBoundary );

      viscousOperator.apply( tmp_, dst, level, flag, updateType );
      projectNormal_.project( dst, level, FreeslipBoundary );
   }

   void computeInverseDiagonalOperatorValues() { viscousOperator.computeInverseDiagonalOperatorValues(); }

   P2VectorFunction< real_t > tmp_;

   ViscousOperator_T        viscousOperator;
   P2ProjectNormalOperator& projectNormal_;
};
} // namespace hyteg