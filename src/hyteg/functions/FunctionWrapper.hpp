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
#include "hyteg/functions/GenericFunction.hpp"

namespace hyteg {

using walberla::uint_t;

/// Base class for all function classes in HyTeG
template < typename func_t >
class FunctionWrapper final : public GenericFunction< typename FunctionTrait< func_t >::ValueType >
{
 public:
   typedef typename FunctionTrait< func_t >::ValueType value_t;

   /// No need for this one, if we do not implement a setter method for wrappedFunc_;
   FunctionWrapper() = delete;

   /// Constructor that constructs the function which the class wraps itself around
   FunctionWrapper( const std::string&                         name,
                    const std::shared_ptr< PrimitiveStorage >& storage,
                    size_t                                     minLevel,
                    size_t                                     maxLevel )
   {
      wrappedFunc_ = std::make_unique< func_t >( name, storage, minLevel, maxLevel );
   };

   ~FunctionWrapper()
   {
     WALBERLA_LOG_INFO_ON_ROOT( "Destructing '" << this->getFunctionName() << "'" );
   }

   /// provide access to wrapped function
   /// @{
   func_t& unwrap() { return *wrappedFunc_; }

   const func_t& unwrap() const { return *wrappedFunc_; }
   /// @}

   uint_t getDimension() const { return wrappedFunc_->getDimension(); };

   const std::string& getFunctionName() const { return wrappedFunc_->getFunctionName(); };

   std::shared_ptr< PrimitiveStorage > getStorage() const { return wrappedFunc_->getStorage(); }

   void multElementwise( const std::vector< std::reference_wrapper< const GenericFunction< value_t > > >& functions,
                         uint_t                                                                           level,
                         DoFType                                                                          flag = All ) const
   {
      std::vector< std::reference_wrapper< const func_t > > realFuncs;
      for ( const GenericFunction< value_t >& func : functions )
      {
         realFuncs.push_back( func.template unwrap< func_t >() );
      }
      wrappedFunc_->multElementwise( realFuncs, level, flag );
   };

   void interpolate( value_t constant, uint_t level, DoFType flag = All ) const
   {
      wrappedFunc_->interpolate( constant, level, flag );
   };

   void interpolate( const std::function< value_t( const hyteg::Point3D& ) >& expr, uint_t level, DoFType flag = All ) const
   {
      wrappedFunc_->interpolate( expr, level, flag );
   };

   void interpolate( const std::vector< std::function< value_t( const hyteg::Point3D& ) > >& expressions,
                     uint_t                                                                  level,
                     DoFType                                                                 flag = All ) const
   {
      wrappedFunc_->interpolate( expressions, level, flag );
   };

   value_t dotGlobal( const GenericFunction< value_t >& secondOp, const uint_t level, const DoFType flag = All ) const
   {
      const func_t& aux = secondOp.template unwrap< func_t >();
      wrappedFunc_->dotGlobal( aux, level, flag );
   };

   value_t dotLocal( const GenericFunction< value_t >& secondOp, uint_t level, DoFType flag = All ) const
   {
      wrappedFunc_->dotLocal( secondOp.template unwrap< func_t >(), level, flag );
   };

   void enableTiming( const std::shared_ptr< walberla::WcTimingTree >& timingTree ) { wrappedFunc_->enableTiming( timingTree ); };

   void setBoundaryCondition( BoundaryCondition bc ) { wrappedFunc_->setBoundaryCondition( bc ); };

   BoundaryCondition getBoundaryCondition() const { return wrappedFunc_->getBoundaryCondition(); };

   void add( const value_t scalar, uint_t level, DoFType flag = All ) const { wrappedFunc_->add( scalar, level, flag ); };

   void add( const std::vector< value_t >                                                     scalars,
             const std::vector< std::reference_wrapper< const GenericFunction< value_t > > >& functions,
             uint_t                                                                           level,
             DoFType                                                                          flag = All ) const
   {
      std::vector< std::reference_wrapper< const func_t > > realFuncs;
      for ( const GenericFunction< value_t >& func : functions )
      {
         realFuncs.push_back( func.template unwrap< func_t >() );
      }
      wrappedFunc_->add( scalars, realFuncs, level, flag );
   };

   void assign( const std::vector< value_t >                                                     scalars,
                const std::vector< std::reference_wrapper< const GenericFunction< value_t > > >& functions,
                uint_t                                                                           level,
                DoFType                                                                          flag = All ) const
   {
      std::vector< std::reference_wrapper< const func_t > > realFuncs;
      for ( const GenericFunction< value_t >& func : functions )
      {
         realFuncs.push_back( func.template unwrap< func_t >() );
      }
      wrappedFunc_->assign( scalars, realFuncs, level, flag );
   };

   void swap( const GenericFunction< value_t >& other, const uint_t& level, const DoFType& flag = All ) const
   {
      wrappedFunc_->swap( other.template unwrap< func_t >(), level, flag );
   };

   void copyFrom( const GenericFunction< value_t >&              other,
                  const uint_t&                                  level,
                  const std::map< PrimitiveID::IDType, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID::IDType, uint_t >& otherPrimitiveIDsToRank ) const
   {
      wrappedFunc_->copyFrom( other.template unwrap< func_t >(), level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
   };

 private:
   std::unique_ptr< func_t > wrappedFunc_;
};

} // namespace hyteg
