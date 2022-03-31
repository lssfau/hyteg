/*
 * Copyright (c) 2017-2020 Marcus Mohr.
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

#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/operators/VectorToVectorOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2SurrogateOperator.hpp"
#include "hyteg/p2functionspace/P2VariableOperator.hpp"
#include "hyteg/solvers/Smoothables.hpp"

namespace hyteg {

using walberla::real_t;

template < typename ValueType, template < typename > class VecFuncKind, class SubOpType >
class VectorLaplaceOperator : public VectorToVectorOperator< ValueType, VecFuncKind, VecFuncKind >,
                              public WeightedJacobiSmoothable< VecFuncKind< ValueType > >,
                              public GSSmoothable< VecFuncKind< ValueType > >,
                              public GSBackwardsSmoothable< VecFuncKind< ValueType > >,
                              public SORSmoothable< VecFuncKind< ValueType > >,
                              public SORBackwardsSmoothable< VecFuncKind< ValueType > >,
                              public OperatorWithInverseDiagonal< VecFuncKind< ValueType > >
{
 public:
   template < typename... SpecialCtorArgs >
   VectorLaplaceOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                          size_t                                     minLevel,
                          size_t                                     maxLevel,
                          SpecialCtorArgs... extraArgs )
   : VectorToVectorOperator< ValueType, VecFuncKind, VecFuncKind >( storage, minLevel, maxLevel )
   {
      std::shared_ptr< SubOpType > zero( nullptr );
      std::shared_ptr< SubOpType > lapl = std::make_shared< SubOpType >( storage, minLevel, maxLevel, extraArgs... );

      if ( this->dim_ == 3 )
      {
         this->subOper_[0][0] = lapl;
         this->subOper_[0][1] = zero;
         this->subOper_[0][2] = zero;

         this->subOper_[1][0] = zero;
         this->subOper_[1][1] = lapl;
         this->subOper_[1][2] = zero;

         this->subOper_[2][0] = zero;
         this->subOper_[2][1] = zero;
         this->subOper_[2][2] = lapl;
      }
      else
      {
         this->subOper_[0][0] = lapl;
         this->subOper_[0][1] = zero;

         this->subOper_[1][0] = zero;
         this->subOper_[1][1] = lapl;
      }
   }

   void smooth_jac( const VecFuncKind< ValueType >& dst,
                    const VecFuncKind< ValueType >& rhs,
                    const VecFuncKind< ValueType >& src,
                    real_t                          relax,
                    size_t                          level,
                    DoFType                         flag ) const final;

   void smooth_gs( const VecFuncKind< ValueType >& dst,
                   const VecFuncKind< ValueType >& rhs,
                   size_t                          level,
                   DoFType                         flag ) const final;

   void smooth_gs_backwards( const VecFuncKind< ValueType >& dst,
                             const VecFuncKind< ValueType >& rhs,
                             size_t                          level,
                             DoFType                         flag ) const final;

   void smooth_sor( const VecFuncKind< ValueType >& dst,
                    const VecFuncKind< ValueType >& rhs,
                    real_t                          relax,
                    size_t                          level,
                    DoFType                         flag ) const final;

   void smooth_sor_backwards( const VecFuncKind< ValueType >& dst,
                              const VecFuncKind< ValueType >& rhs,
                              real_t                          relax,
                              size_t                          level,
                              DoFType                         flag ) const final;

   std::shared_ptr< VecFuncKind< ValueType > > getInverseDiagonalValues() const final
   {
      return this->extractInverseDiagonal();
   }

   void computeInverseDiagonalOperatorValues() final
   {
      this->VectorToVectorOperator< ValueType, VecFuncKind, VecFuncKind >::computeInverseDiagonalOperatorValues();
   }

};

// ------------------------
//  stencil-based versions
// ------------------------
typedef VectorLaplaceOperator< real_t, P1VectorFunction, P1ConstantLaplaceOperator > P1ConstantVectorLaplaceOperator;
typedef VectorLaplaceOperator< real_t, P2VectorFunction, P2ConstantLaplaceOperator > P2ConstantVectorLaplaceOperator;

// ----------------------
//  elementwise versions
// ----------------------
typedef VectorLaplaceOperator< real_t, P1VectorFunction, P1ElementwiseLaplaceOperator > P1ElementwiseVectorLaplaceOperator;
typedef VectorLaplaceOperator< real_t, P1VectorFunction, P1ElementwiseBlendingLaplaceOperator >
    P1ElementwiseBlendingVectorLaplaceOperator;

typedef VectorLaplaceOperator< real_t, P2VectorFunction, P2ElementwiseLaplaceOperator > P2ElementwiseVectorLaplaceOperator;
typedef VectorLaplaceOperator< real_t, P2VectorFunction, P2ElementwiseBlendingLaplaceOperator >
    P2ElementwiseBlendingVectorLaplaceOperator;

// ------------------
//  special versions
// ------------------
typedef VectorLaplaceOperator< real_t, P2VectorFunction, P2BlendingLaplaceOperator >  P2BlendingVectorLaplaceOperator;
typedef VectorLaplaceOperator< real_t, P2VectorFunction, P2SurrogateLaplaceOperator > P2SurrogateVectorLaplaceOperator;

} // namespace hyteg
