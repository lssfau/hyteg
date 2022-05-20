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

#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/functions/FunctionWrapper.hpp"

namespace hyteg {

using walberla::uint_t;

/// Base class for Block Functions
///
/// The implicit assumptions here are that
/// - all sub-functions have the same value type
/// - all sub-functions work on the same primitive storage
/// - all sub-functions live on the same range of levels
template < typename value_t >
class BlockFunction
{
 public:
   typedef value_t valueType;

   template < typename VType >
   using FunctionType = BlockFunction< VType >;

   BlockFunction( const std::string name )
   : functionName_( name )
   {}

   /// @name Query Functions
   /// Methods for questioning object for certain properties
   /// @{
   const std::string& getFunctionName() const { return functionName_; }

   std::shared_ptr< PrimitiveStorage > getStorage() const { return subFunc_[0]->getStorage(); }

   /// A block function in itself has no dimension (only its sub-functions do)
   uint_t getDimension() const { return 0; }

   /// Query the number of "blocks", i.e. number of sub-functions
   uint_t getNumberOfBlocks() const { return subFunc_.size(); }
   /// @}

   /// @name Component access
   /// Methods for component function access
   /// @{
   GenericFunction< value_t >& getSubFunction( uint_t idx )
   {
      WALBERLA_ASSERT_LESS( idx, subFunc_.size() );
      return *( subFunc_[idx] );
   }

   const GenericFunction< value_t >& getSubFunction( uint_t idx ) const
   {
      WALBERLA_ASSERT_LESS( idx, subFunc_.size() );
      return *( subFunc_[idx] );
   }

   GenericFunction< value_t >& operator[]( uint_t idx ) { return getSubFunction( idx ); }

   const GenericFunction< value_t >& operator[]( uint_t idx ) const { return getSubFunction( idx ); }
   /// @}

   void multElementwise( const std::vector< std::reference_wrapper< const BlockFunction< value_t > > >& functions,
                         uint_t                                                                         level,
                         DoFType                                                                        flag = All ) const
   {
      WALBERLA_DEBUG_SECTION()
      {
         for ( auto func : functions )
         {
            WALBERLA_ASSERT_EQUAL( subFunc_.size(), func.get().getNumberOfBlocks() );
         }
      }

      for ( uint_t k = 0; k < subFunc_.size(); k++ )
      {
         subFunc_[k]->multElementwise( filter( k, functions ), level, flag );
      }
   }

   void interpolate( const std::function< value_t( const hyteg::Point3D& ) >& expr, size_t level, DoFType flag = All ) const
   {
      for ( uint_t k = 0; k < subFunc_.size(); k++ )
      {
         subFunc_[k]->interpolate( expr, level, flag );
      }
   }

   void interpolate( value_t constant, size_t level, DoFType flag = All ) const
   {
      for ( uint_t k = 0; k < subFunc_.size(); k++ )
      {
         subFunc_[k]->interpolate( constant, level, flag );
      }
   }

   void interpolate( const std::vector< std::function< value_t( const hyteg::Point3D& ) > >& expr,
                     size_t                                                                  level,
                     DoFType                                                                 flag = All ) const
   {
      WALBERLA_ASSERT_GREATER_EQUAL( expr.size(), subFunc_.size() );
      WALBERLA_DEBUG_SECTION()
      {
         if ( expr.size() > subFunc_.size() )
         {
            WALBERLA_LOG_WARNING_ON_ROOT( "BlockFunction::interpolate(): Ignoring excess std::functions!" );
         }
      }
      for ( uint_t k = 0; k < subFunc_.size(); ++k )
      {
         subFunc_[k]->interpolate( expr[k], level, flag );
      }
   }

   void add( value_t scalar, size_t level, DoFType flag = All ) const
   {
      for ( uint_t k = 0; k < subFunc_.size(); ++k )
      {
         subFunc_[k]->add( scalar, level, flag );
      }
   }

   void add( const std::vector< value_t >                                                   scalars,
             const std::vector< std::reference_wrapper< const BlockFunction< value_t > > >& functions,
             size_t                                                                         level,
             DoFType                                                                        flag = All ) const
   {
      WALBERLA_DEBUG_SECTION()
      {
         for ( auto func : functions )
         {
            WALBERLA_ASSERT_EQUAL( subFunc_.size(), func.get().getNumberOfBlocks() );
         }
      }

      for ( uint_t k = 0; k < subFunc_.size(); ++k )
      {
         subFunc_[k]->add( scalars, filter( k, functions ), level, flag );
      }
   }

   void enableTiming( const std::shared_ptr< walberla::WcTimingTree >& timingTree )
   {
      for ( uint_t k = 0; k < subFunc_.size(); ++k )
      {
         subFunc_[k]->enableTiming( timingTree );
      }
   }

   void swap( const BlockFunction< value_t >& other, const uint_t& level, const DoFType& flag = All ) const
   {
      WALBERLA_ASSERT_EQUAL( subFunc_.size(), other.subFunc_.size() );
      for ( uint_t k = 0; k < subFunc_.size(); ++k )
      {
         subFunc_[k]->swap( other[k], level, flag );
      }
   }

   void assign( const std::vector< value_t >&                                                  scalars,
                const std::vector< std::reference_wrapper< const BlockFunction< value_t > > >& functions,
                size_t                                                                         level,
                DoFType                                                                        flag = All ) const
   {
      WALBERLA_DEBUG_SECTION()
      {
         for ( auto func : functions )
         {
            WALBERLA_ASSERT_EQUAL( subFunc_.size(), func.get().getNumberOfBlocks() );
         }
      }

      for ( uint_t k = 0; k < subFunc_.size(); ++k )
      {
         subFunc_[k]->assign( scalars, filter( k, functions ), level, flag );
      }
   }

   value_t dotLocal( const BlockFunction< value_t >& rhs, const uint_t level, const DoFType flag = All ) const
   {
      value_t sum = value_t( 0 );
      for ( uint_t k = 0; k < subFunc_.size(); ++k )
      {
         sum += subFunc_[k]->dotLocal( rhs[k], level, flag );
      }
      return sum;
   }

   value_t dotGlobal( const BlockFunction< value_t >& rhs, const uint_t level, const DoFType flag = All ) const
   {
      auto sum = dotLocal( rhs, level, flag );
      walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
      return sum;
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
   void copyFrom( const BlockFunction< value_t >&                other,
                  const uint_t&                                  level,
                  const std::map< PrimitiveID::IDType, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID::IDType, uint_t >& otherPrimitiveIDsToRank ) const
   {
      for ( uint_t k = 0; k < subFunc_.size(); ++k )
      {
         subFunc_[k]->copyFrom( other[k], level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
      }
   }

   uint_t getNumberOfLocalDoFs( uint_t level ) const
   {
      uint_t nDoFs = 0;
      for ( uint_t k = 0; k < subFunc_.size(); ++k )
      {
         nDoFs += subFunc_[k]->getNumberOfLocalDoFs( level );
      }
      return nDoFs;
   }

   uint_t getNumberOfGlobalDoFs( uint_t          level,
                                 const MPI_Comm& communicator = walberla::mpi::MPIManager::instance()->comm(),
                                 const bool&     onRootOnly   = false ) const
   {
      uint_t nDoFs = 0;
      for ( uint_t k = 0; k < subFunc_.size(); ++k )
      {
         nDoFs += subFunc_[k]->getNumberOfGlobalDoFs( level, communicator, onRootOnly );
      }
      return nDoFs;
   }

   void enumerate( uint_t level ) const
   {
      uint_t counterDoFs = getNumberOfLocalDoFs( level );

      std::vector< uint_t > doFsPerRank = walberla::mpi::allGather( counterDoFs );

      valueType offset = 0;

      for ( uint_t i = 0; i < uint_c( walberla::MPIManager::instance()->rank() ); ++i )
      {
         offset += static_cast< valueType >( doFsPerRank[i] );
      }

      for ( uint_t k = 0; k < subFunc_.size(); k++ )
      {
         subFunc_[k]->enumerate( level, offset );
      }
   }

   template < typename OtherType >
   void copyBoundaryConditionFromFunction( const BlockFunction< OtherType >& other )
   {
      WALBERLA_ASSERT_EQUAL( subFunc_.size(), other.getNumberOfBlocks() );
      for ( uint_t k = 0; k < subFunc_.size(); k++ )
      {
         subFunc_[k]->setBoundaryCondition( other.getSubFunction( k ).getBoundaryCondition() );
      }
   }

   /// conversion to/from linear algebra representation
   /// @{
   void toVector( const BlockFunction< idx_t >&         numerator,
                  const std::shared_ptr< VectorProxy >& vec,
                  uint_t                                level,
                  DoFType                               flag ) const
   {
      for ( uint_t k = 0; k < subFunc_.size(); ++k )
      {
         subFunc_[k]->toVector( numerator[k], vec, level, flag );
      }
   }

   void fromVector( const BlockFunction< idx_t >&         numerator,
                    const std::shared_ptr< VectorProxy >& vec,
                    uint_t                                level,
                    DoFType                               flag ) const
   {
      for ( uint_t k = 0; k < subFunc_.size(); ++k )
      {
         subFunc_[k]->fromVector( numerator[k], vec, level, flag );
      }
   };
   /// @}

 protected:
   const std::string                                            functionName_;
   std::vector< std::shared_ptr< GenericFunction< value_t > > > subFunc_;

   // extract a vector of generic sub-functions from a vector of block functions
   std::vector< std::reference_wrapper< const GenericFunction< value_t > > >
       filter( uint_t idx, const std::vector< std::reference_wrapper< const BlockFunction< value_t > > >& functions ) const
   {
      std::vector< std::reference_wrapper< const GenericFunction< value_t > > > subFuncVec;
      std::transform( functions.begin(), functions.end(), std::back_inserter( subFuncVec ), [idx]( auto& function ) {
         return std::cref( function.get().getSubFunction( idx ) );
      } );
      return subFuncVec;
   }
};

} // namespace hyteg
