/*
 * Copyright (c) 2017-2021 Marcus Mohr.
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
#include "hyteg/operators/VectorLaplaceOperator.hpp"

namespace hyteg {

using walberla::real_t;

template < typename ValueType, template < typename > class VecFuncKind, class SubOpType >
void VectorLaplaceOperator< ValueType, VecFuncKind, SubOpType >::smooth_jac( const VecFuncKind< ValueType >& dst,
                                                                             const VecFuncKind< ValueType >& rhs,
                                                                             const VecFuncKind< ValueType >& src,
                                                                             real_t                          relax,
                                                                             size_t                          level,
                                                                             DoFType                         flag ) const
{
   for ( uint_t k = 0; k < this->dim_; ++k )
   {
      if ( const auto subOp =
               std::dynamic_pointer_cast< WeightedJacobiSmoothable< typename SubOpType::srcType > >( this->subOper_[k][k] ) )
      {
         subOp->smooth_jac( dst[k], rhs[k], src[k], relax, level, flag );
      }
      else
      {
         throw std::runtime_error(
             "Jacobi smoothing of a VectorLaplaceOperator requires its diagonal blocks to have the WeightedJacobiSmoothable interface." );
      }
   }
}

template < typename ValueType, template < typename > class VecFuncKind, class SubOpType >
void VectorLaplaceOperator< ValueType, VecFuncKind, SubOpType >::smooth_gs( const VecFuncKind< ValueType >& dst,
                                                                            const VecFuncKind< ValueType >& rhs,
                                                                            size_t                          level,
                                                                            DoFType                         flag ) const
{
   for ( uint_t k = 0; k < this->dim_; ++k )
   {
      if ( const auto* subOp = dynamic_cast< const GSSmoothable< typename SubOpType::srcType >* >( this->subOper_[k][k].get() ) )
      {
         subOp->smooth_gs( dst[k], rhs[k], level, flag );
      }
      else
      {
         throw std::runtime_error(
             "Gauss-Seidel smoothing of a VectorLaplaceOperator requires its diagonal blocks to have the GSSmoothable interface." );
      }
   }
}

template < typename ValueType, template < typename > class VecFuncKind, class SubOpType >
void VectorLaplaceOperator< ValueType, VecFuncKind, SubOpType >::smooth_gs_backwards( const VecFuncKind< ValueType >& dst,
                                                                                      const VecFuncKind< ValueType >& rhs,
                                                                                      size_t                          level,
                                                                                      DoFType                         flag ) const
{
   for ( uint_t k = 0; k < this->dim_; ++k )
   {
      if ( const auto* subOp =
               dynamic_cast< const GSBackwardsSmoothable< typename SubOpType::srcType >* >( this->subOper_[k][k].get() ) )
      {
         subOp->smooth_gs_backwards( dst[k], rhs[k], level, flag );
      }
      else
      {
         throw std::runtime_error(
             "Backwards Gauss-Seidel smoothing of a VectorLaplaceOperator requires its diagonal blocks to have the GSBackwardsSmoothable interface." );
      }
   }
}

template < typename ValueType, template < typename > class VecFuncKind, class SubOpType >
void VectorLaplaceOperator< ValueType, VecFuncKind, SubOpType >::smooth_sor( const VecFuncKind< ValueType >& dst,
                                                                             const VecFuncKind< ValueType >& rhs,
                                                                             real_t                          relax,
                                                                             size_t                          level,
                                                                             DoFType                         flag ) const
{
   for ( uint_t k = 0; k < this->dim_; ++k )
   {
      if ( const auto* subOp = dynamic_cast< const SORSmoothable< typename SubOpType::srcType >* >( this->subOper_[k][k].get() ) )
      {
         subOp->smooth_sor( dst[k], rhs[k], relax, level, flag );
      }
      else
      {
         throw std::runtime_error(
             "SOR smoothing of a VectorLaplaceOperator requires its diagonal blocks to have the SORSmoothable interface." );
      }
   }
}

template < typename ValueType, template < typename > class VecFuncKind, class SubOpType >
void VectorLaplaceOperator< ValueType, VecFuncKind, SubOpType >::smooth_sor_backwards( const VecFuncKind< ValueType >& dst,
                                                                                       const VecFuncKind< ValueType >& rhs,
                                                                                       real_t                          relax,
                                                                                       size_t                          level,
                                                                                       DoFType flag ) const
{
   for ( uint_t k = 0; k < this->dim_; ++k )
   {
      if ( const auto* subOp =
               dynamic_cast< const SORBackwardsSmoothable< typename SubOpType::srcType >* >( this->subOper_[k][k].get() ) )
      {
         subOp->smooth_sor_backwards( dst[k], rhs[k], relax, level, flag );
      }
      else
      {
         throw std::runtime_error(
             "Backward SOR smoothing of a VectorLaplaceOperator requires its diagonal blocks to have the SORBackwardsSmoothable interface." );
      }
   }
}

// =========================
//  Explicit Instantiations
// =========================

// P1ConstantVectorLaplaceOperator
template class VectorLaplaceOperator< real_t, P1VectorFunction, P1ConstantLaplaceOperator >;

// P1ElementwiseVectorLaplaceOperator
template class VectorLaplaceOperator< real_t, P1VectorFunction, P1ElementwiseLaplaceOperator >;

// P1ElementwiseBlendingVectorLaplaceOperator
template class VectorLaplaceOperator< real_t, P1VectorFunction, P1ElementwiseBlendingLaplaceOperator >;

// P2ConstantVectorLaplaceOperator
template class VectorLaplaceOperator< real_t, P2VectorFunction, P2ConstantLaplaceOperator >;

// P2ElementwiseVectorLaplaceOperator
template class VectorLaplaceOperator< real_t, P2VectorFunction, P2ElementwiseLaplaceOperator >;

// P2ElementwiseBlendingVectorLaplaceOperator
template class VectorLaplaceOperator< real_t, P2VectorFunction, P2ElementwiseBlendingLaplaceOperator >;

// P2BlendingVectorLaplaceOperator
template class VectorLaplaceOperator< real_t, P2VectorFunction, P2BlendingLaplaceOperator >;

// P2SurrogateVectorLaplaceOperator
template class VectorLaplaceOperator< real_t, P2VectorFunction, P2SurrogateLaplaceOperator >;

} // namespace hyteg
