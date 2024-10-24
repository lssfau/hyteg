/*
 * Copyright (c) 2021 Marcus Mohr.
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

#include "core/DataTypes.h"

#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/operators/GenericOperator.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"

namespace hyteg {

using walberla::uint_t;

template < typename oper_t >
class OperatorWrapper final
: public GenericOperator< typename FunctionTrait< typename oper_t::srcType >::ValueType >,
  public GSSmoothable< GenericFunction< typename FunctionTrait< typename oper_t::srcType >::ValueType > >
{
 public:
   typedef typename FunctionTrait< typename oper_t::srcType >::ValueType value_t;

   /// No need for this one, if we do not implement a setter method for wrappedFunc_;
   OperatorWrapper() = delete;

   /// Constructor that constructs the operator which the class wraps itself around
   template < typename... ConstructorArguments >
   OperatorWrapper( const std::shared_ptr< PrimitiveStorage >& storage,
                    uint_t                                     minLevel,
                    uint_t                                     maxLevel,
                    ConstructorArguments... args )
   {
      wrappedOper_ = std::make_unique< oper_t >( storage, minLevel, maxLevel, args... );
   }

   ~OperatorWrapper() {}

   /// provide access to wrapped operator
   /// @{
   oper_t& unwrap() { return *wrappedOper_; }

   const oper_t& unwrap() const { return *wrappedOper_; }
   /// @}

   void apply( const GenericFunction< value_t >& src,
               const GenericFunction< value_t >& dst,
               size_t                            level,
               DoFType                           flag,
               UpdateType                        updateType = Replace ) const
   {
      wrappedOper_->apply( src.template unwrap< typename oper_t::srcType >(),
                           dst.template unwrap< typename oper_t::dstType >(),
                           level,
                           flag,
                           updateType );
   };

   virtual void gemv( const value_t&                    alpha,
                      const GenericFunction< value_t >& src,
                      const value_t&                    beta,
                      const GenericFunction< value_t >& dst,
                      size_t                            level,
                      DoFType                           flag ) const
   {
      wrappedOper_->gemv( alpha,
                          src.template unwrap< typename oper_t::srcType >(),
                          beta,
                          dst.template unwrap< typename oper_t::dstType >(),
                          level,
                          flag );
   }

   void smooth_gs( const GenericFunction< value_t >& dst,
                   const GenericFunction< value_t >& rhs,
                   size_t                            level,
                   DoFType                           flag ) const override
   {
      // the "real" types of our operator
      using srcType = typename oper_t::srcType;
      using dstType = typename oper_t::dstType;

      // if both operator types are not the same we will violate the GSSmoothable interface
      // we throw an exception in this case, since GS is not defined properly in this case
      if constexpr ( std::is_same< srcType, dstType >::value )
      {
         if ( const auto* op_gs = dynamic_cast< const GSSmoothable< srcType >* >( wrappedOper_.get() ) )
         {
            const auto& dst_unwrapped = dst.template unwrap< typename oper_t::srcType >();
            const auto& rhs_unwrapped = rhs.template unwrap< typename oper_t::dstType >();
            op_gs->smooth_gs( dst_unwrapped, rhs_unwrapped, level, flag );
         }
         else
         {
            throw std::runtime_error( "The Gauss-Seidel Operator requires the GSSmoothable interface." );
         }
      }
      else
      {
         throw std::runtime_error( "For GaussSeidel src and dst functions must coincide" );
      }
   }

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const GenericFunction< idx_t >&             src,
                  const GenericFunction< idx_t >&             dst,
                  size_t                                      level,
                  DoFType                                     flag ) const
   {
      wrappedOper_->toMatrix( mat,
                              src.template unwrap< typename oper_t::srcType::template FunctionType< idx_t > >(),
                              dst.template unwrap< typename oper_t::dstType::template FunctionType< idx_t > >(),
                              level,
                              flag );
   };

 private:
   std::unique_ptr< oper_t > wrappedOper_;
};

} // namespace hyteg
