/*
 * Copyright (c) 2020 Marcus Mohr.
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

#include "hyteg/functions/Function.hpp"
#include "hyteg/functions/GenericFunction.hpp"

namespace hyteg {

using walberla::uint_t;

/// Base class for VectorFunctions
///
/// This is the base class for all function classes in HyTeG representing vector fileds
/// which have scalar component functions, i.e. where the VectorClass is basically a
/// Container for Scalar Functions (CSF)
template < typename ValueType >
class CSFVectorFunction : public GenericFunction
{
 public:
   CSFVectorFunction( const std::string name )
   : functionName_( name )
   {}

   const std::string& getFunctionName() const { return functionName_; }

   std::shared_ptr< PrimitiveStorage > getStorage() const { return compFunc_[0]->getStorage(); }

   /// \note Dimension of VectorFunction is decoupled from storage now
   uint_t getDimension() const { return compFunc_.size(); }

   /// @name Component access
   /// Methods for component function access
   /// @{
   Function< ValueType >& component( uint_t idx )
   {
      WALBERLA_ASSERT_LESS( idx, compFunc_.size() );
      return *( compFunc_[idx] );
   }

   const Function< ValueType >& component( uint_t idx ) const
   {
      WALBERLA_ASSERT_LESS( idx, compFunc_.size() );
      return *( compFunc_[idx] );
   }

   Function< ValueType >& operator[]( uint_t idx ) { return component( idx ); }

   const Function< ValueType >& operator[]( uint_t idx ) const { return component( idx ); }
   /// @}

   void multElementwise( const std::vector< std::reference_wrapper< const CSFVectorFunction< ValueType > > >& functions,
                         uint_t                                                                               level,
                         DoFType                                                                              flag = All ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); k++ )
      {
         compFunc_[k]->multElementwise( vectorFunctionTools::filter( 0, functions ), level, flag );
      }
   }

   void interpolate( const std::function< ValueType( const hyteg::Point3D& ) >& expr, size_t level, DoFType flag = All ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); k++ )
      {
         compFunc_[k]->interpolate( expr, level, flag );
      }
   }

   void interpolate( const ValueType& constant, size_t level, DoFType flag = All ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); k++ )
      {
         compFunc_[k]->interpolate( constant, level, flag );
      }
   }

   void interpolate( const std::vector< std::function< ValueType( const hyteg::Point3D& ) > >& expr,
                     size_t                                                                    level,
                     DoFType                                                                   flag = All ) const
   {
      WALBERLA_ASSERT_GREATER( expr.size(), 0 );
      WALBERLA_ASSERT_EQUAL( expr.size(), compFunc_.size() );

      for ( uint_t k = 0; k < expr.size(); ++k )
      {
         compFunc_[k].interpolate( expr[k], level, flag );
      }
   }

   void add( real_t scalar, size_t level, DoFType flag = All ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k].add( scalar, level, flag );
      }
   }

   void add( const std::vector< walberla::real_t >                                                scalars,
             const std::vector< std::reference_wrapper< const CSFVectorFunction< ValueType > > >& functions,
             size_t                                                                               level,
             DoFType                                                                              flag = All ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k].add( scalars, vectorFunctionTools::filter( k, functions ), level, flag );
      }
   }

   void enableTiming( const std::shared_ptr< walberla::WcTimingTree >& timingTree )
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->enableTiming( timingTree );
      }
   }

   template < typename OtherFunctionValueType >
   void copyBoundaryConditionFromFunction( const CSFVectorFunction< OtherFunctionValueType >& other )
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->setBoundaryCondition( vectorFunctionTools::filter( k, other ).getBoundaryCondition() );
      }
   }

   void swap( const CSFVectorFunction< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->swap( other.u, level, flag );
      }
   }

   void assign( const std::vector< walberla::real_t >                                                scalars,
                const std::vector< std::reference_wrapper< const CSFVectorFunction< ValueType > > >& functions,
                size_t                                                                               level,
                DoFType                                                                              flag = All ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->assign( scalars, vectorFunctionTools::filter( k, functions ), level, flag );
      }
   }

   walberla::real_t dotLocal( const CSFVectorFunction< ValueType >& rhs, const uint_t level, const DoFType flag = All ) const
   {
      ValueType sum = ValueType( 0 );
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         sum += compFunc_[k]->dotLocal( *rhs[k], level, flag );
      }
      return sum;
   }

   walberla::real_t dotGlobal( const CSFVectorFunction< ValueType >& rhs, const uint_t level, const DoFType flag = All ) const
   {
      auto sum = dotLocal( rhs, level, flag );
      walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
      return sum;
   }

   ValueType getMaxComponentMagnitude( uint_t level, DoFType flag, bool mpiReduce = true ) const
   {
      std::vector< ValueType > values;

      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         values.push_back( compFunc_[k]->getMaxMagnitude( level, flag, mpiReduce ) );
      }

      return *std::max_element( values.begin(), values.end() );
   }

   /// \brief Copies all values function data from other to this.
   ///
   /// This method can be used safely if the other function is located on a different PrimitiveStorage.
   /// This method also works, if the storages are distributed differently.
   ///
   /// \param other another function
   /// \param level the refinement level
   /// \param localPrimitiveIDsToRank Map that contains as keys all primitive IDs of all primitives that are local regarding the
   ///                                storage of this function, and as values the MPI ranks of the processes that own these
   ///                                primitives regarding the storage of the other function
   /// \param otherPrimitiveIDsToRank Map that contains as keys all primitive IDs of all primitives that are local regarding the
   ///                                storage of the other function, and as values the MPI ranks of the processes that own these
   ///                                primitives regarding the storage this function lives on.
   ///
   void copyFrom( const CSFVectorFunction< ValueType >&          other,
                  const uint_t&                                  level,
                  const std::map< PrimitiveID::IDType, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID::IDType, uint_t >& otherPrimitiveIDsToRank ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->copyFrom(
             vectorFunctionTools::filter( k, other ), level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
      }
   }

 protected:
   const std::string                                       functionName_;
   std::vector< std::shared_ptr< Function< ValueType > > > compFunc_;
};

} // namespace hyteg
