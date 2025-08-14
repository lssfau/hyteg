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

#include "hyteg/p2functionspace/P2RotationOperator.hpp"
#include "hyteg_operators/operators/full_stokes/P2VectorElementwiseFullStokesP0ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/full_stokes/P2VectorElementwiseFullStokesP1ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

namespace hyteg {

template < typename ViscosityFunction_T, typename ViscousOperator_T >
class P2ABlockViscousRotationFSTemplate : public Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >,
                                          public OperatorWithInverseDiagonal< P2VectorFunction< real_t > >
{
 public:
   P2ABlockViscousRotationFSTemplate( const std::shared_ptr< PrimitiveStorage >&    storage,
                                      uint_t                                        minLevel,
                                      uint_t                                        maxLevel,
                                      const ViscosityFunction_T&                    mu,
                                      P2RotationOperator&                           rotationOperator,
                                      BoundaryCondition                             bcVelocity,
                                      std::shared_ptr< P2VectorFunction< real_t > > tmp         = nullptr,
                                      std::shared_ptr< P2VectorFunction< real_t > > tmpdst      = nullptr,
                                      std::shared_ptr< P2VectorFunction< idx_t > >  tmpAssembly = nullptr )
   : Operator< P2VectorFunction< real_t >, P2VectorFunction< real_t > >( storage, minLevel, maxLevel )
   , viscousOperator_( storage, minLevel, maxLevel, mu )
   , rotationOperator_( rotationOperator )
   {
      if ( tmp == nullptr )
      {
         tmp_ = std::make_shared< P2VectorFunction< real_t > >(
             "tmp__P2ABlockViscousRotationFSTemplate", storage, minLevel, maxLevel, bcVelocity );
      }
      else
      {
         tmp_ = tmp;
      }

      if ( tmpdst == nullptr )
      {
         tmpdst_ = std::make_shared< P2VectorFunction< real_t > >(
             "tmpdst__P2ABlockViscousRotationFSTemplate", storage, minLevel, maxLevel, bcVelocity );
      }
      else
      {
         tmpdst_ = tmpdst;
      }

      if ( tmpAssembly == nullptr )
      {
         tmpAssembly_ = std::make_shared< P2VectorFunction< idx_t > >(
             "tmpAssembly__P2ABlockViscousRotationFSTemplate", storage, minLevel, maxLevel, bcVelocity );
      }
      else
      {
         tmpAssembly_ = tmpAssembly;
      }
   }

   void apply( const P2VectorFunction< real_t >& src,
               const P2VectorFunction< real_t >& dst,
               const uint_t                      level,
               const DoFType                     flag,
               const UpdateType                  updateType = Replace ) const override
   {
      tmp_->assign( { 1.0 }, { src }, level, All );
      rotationOperator_.rotate( *tmp_, level, FreeslipBoundary, true );

      {
         viscousOperator_.apply( *tmp_, *tmpdst_, level, flag, updateType );
      }

      rotationOperator_.rotate( *tmpdst_, level, FreeslipBoundary );
      dst.assign( { 1.0 }, { *tmpdst_ }, level, flag );
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2VectorFunction< idx_t >&            src,
                  const P2VectorFunction< idx_t >&            dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      auto matProxyOp = mat->createCopy();
      viscousOperator_.toMatrix( matProxyOp, src, dst, level, flag );

      auto matProxyProjectionPost = mat->createCopy();

      rotationOperator_.toMatrix( matProxyProjectionPost, *tmpAssembly_, level, FreeslipBoundary, false );

      std::vector< std::shared_ptr< SparseMatrixProxy > > matrices;
      matrices.push_back( matProxyProjectionPost );
      matrices.push_back( matProxyOp );

      auto matProxyProjectionPre = mat->createCopy();
      rotationOperator_.toMatrix( matProxyProjectionPre, *tmpAssembly_, level, FreeslipBoundary, true );

      matrices.push_back( matProxyProjectionPre );
      mat->createFromMatrixProduct( matrices );
   }

   void computeInverseDiagonalOperatorValues() { viscousOperator_.computeInverseDiagonalOperatorValues(); }

   std::shared_ptr< P2VectorFunction< real_t > > getInverseDiagonalValues() const
   {
      return viscousOperator_.getInverseDiagonalValues();
   }

   std::shared_ptr< P2VectorFunction< real_t > > tmp_;
   std::shared_ptr< P2VectorFunction< real_t > > tmpdst_;
   std::shared_ptr< P2VectorFunction< idx_t > >  tmpAssembly_;

   ViscousOperator_T   viscousOperator_;
   P2RotationOperator& rotationOperator_;
};

using P2ABlockOperatorWithRotation =
    P2ABlockViscousRotationFSTemplate< P2Function< real_t >, operatorgeneration::P2ViscousBlockFullIcosahedralShellMapOperator >;

using P2ABlockP1ViscousOperatorWithRotation =
    P2ABlockViscousRotationFSTemplate< P1Function< real_t >,
                                       operatorgeneration::P2VectorElementwiseFullStokesP1ViscosityIcosahedralShellMap >;

using P2ABlockP0ViscousOperatorWithRotation =
    P2ABlockViscousRotationFSTemplate< P0Function< real_t >,
                                       operatorgeneration::P2VectorElementwiseFullStokesP0ViscosityIcosahedralShellMap >;
} // namespace hyteg