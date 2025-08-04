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
#include "hyteg_operators/operators/epsilon/P2VectorElementwiseEpsilonP1ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/full_stokes/P2VectorElementwiseFullStokesP0ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/full_stokes/P2VectorElementwiseFullStokesP1ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/full_stokes/P2VectorElementwiseFullStokesRotationP0ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

namespace hyteg {

/***************************************************************************************************
NOTE: Here FS denotes FreeSlip, Stokes A block operator is wrapped with a FreeSlip Projection Wrapper
      Changes the linear system from $Ku    = f $
                                  to $PKP^T = Pf$
***************************************************************************************************/

template < typename ViscosityFunction_T, typename ViscousOperator_T >
class P2ABlockViscousProjectionFSTemplate : public Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >,
                                            public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >
{
 public:
   P2ABlockViscousProjectionFSTemplate( const std::shared_ptr< PrimitiveStorage >&    storage,
                                        uint_t                                        minLevel,
                                        uint_t                                        maxLevel,
                                        const ViscosityFunction_T&                    mu,
                                        P2ProjectNormalOperator&                      projectNormal,
                                        BoundaryCondition                             bcVelocity,
                                        std::shared_ptr< P2VectorFunction< real_t > > tmp = nullptr )
   : Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >( storage, minLevel, maxLevel )
   , viscousOperator( storage, minLevel, maxLevel, mu )
   , projectNormal_( projectNormal )
   {
      if ( tmp == nullptr )
      {
         tmp_ = std::make_shared< P2VectorFunction< real_t > >(
             "tmp__P2ViscousIcosahedralShellMapOperatorFS", storage, minLevel, maxLevel, bcVelocity );
      }
      else
      {
         tmp_ = tmp;
      }
   }

   void apply( const P2VectorFunction< real_t >& src,
               const P2VectorFunction< real_t >& dst,
               const uint_t                      level,
               const DoFType                     flag,
               const UpdateType                  updateType = Replace ) const override
   {
      tmp_->assign( { 1 }, { src }, level, All );

      viscousOperator.apply( *tmp_, dst, level, flag, updateType );
      projectNormal_.project( dst, level, FreeslipBoundary );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2VectorFunction< idx_t >&            src,
                  const P2VectorFunction< idx_t >&            dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      viscousOperator.toMatrix( mat, src, dst, level, flag );
   }

   void computeInverseDiagonalOperatorValues() override { viscousOperator.computeInverseDiagonalOperatorValues(); }

   std::shared_ptr< P2VectorFunction< real_t > > getInverseDiagonalValues() const override
   {
      return viscousOperator.getInverseDiagonalValues();
   }

   std::shared_ptr< P2VectorFunction< real_t > > tmp_;

   ViscousOperator_T        viscousOperator;
   P2ProjectNormalOperator& projectNormal_;
};

using P2ABlockOperatorWithProjection =
    P2ABlockViscousProjectionFSTemplate< P2Function< real_t >,
                                         operatorgeneration::P2ViscousBlockFullIcosahedralShellMapOperator >;
using P2ABlockP1ViscousOperatorWithProjection =
    P2ABlockViscousProjectionFSTemplate< P1Function< real_t >,
                                         operatorgeneration::P2VectorElementwiseFullStokesP1ViscosityIcosahedralShellMap >;

using P2ABlockP0ViscousOperatorWithProjection =
    P2ABlockViscousProjectionFSTemplate< P0Function< real_t >,
                                         operatorgeneration::P2VectorElementwiseFullStokesP0ViscosityIcosahedralShellMap >;

} // namespace hyteg