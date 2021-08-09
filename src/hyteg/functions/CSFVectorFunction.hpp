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
#include "hyteg/functions/VectorFunctionTools.hpp"

namespace hyteg {

using walberla::uint_t;

/// Base class for VectorFunctions
///
/// This is the base class for all function classes in HyTeG representing vector fields
/// which have scalar component functions, i.e. where the VectorClass is basically a
/// Container for Scalar Functions (CSF)
template < typename VectorFunctionType >
class CSFVectorFunction
{
 public:
   typedef typename FunctionTrait< VectorFunctionType >::Tag                 Tag;
   typedef typename FunctionTrait< VectorFunctionType >::ValueType           valueType;
   typedef typename FunctionTrait< VectorFunctionType >::VectorComponentType VectorComponentType;

   CSFVectorFunction( const std::string name )
   : functionName_( name )
   {}

   /// @name Query Functions
   /// Methods for questioning object for certain properties
   /// @{
   const std::string& getFunctionName() const { return functionName_; }

   std::shared_ptr< PrimitiveStorage > getStorage() const { return compFunc_[0]->getStorage(); }

   /// \note Dimension of VectorFunction is decoupled from storage now
   uint_t getDimension() const { return compFunc_.size(); }
   /// @}

   /// @name Component access
   /// Methods for component function access
   /// @{
   VectorComponentType& component( uint_t idx )
   {
      WALBERLA_ASSERT_LESS( idx, compFunc_.size() );
      return *( compFunc_[idx] );
   }

   const VectorComponentType& component( uint_t idx ) const
   {
      WALBERLA_ASSERT_LESS( idx, compFunc_.size() );
      return *( compFunc_[idx] );
   }

   VectorComponentType& operator[]( uint_t idx ) { return component( idx ); }

   const VectorComponentType& operator[]( uint_t idx ) const { return component( idx ); }
   /// @}

   void multElementwise( const std::vector< std::reference_wrapper< const VectorFunctionType > >& functions,
                         uint_t                                                                   level,
                         DoFType                                                                  flag = All ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); k++ )
      {
         compFunc_[k]->multElementwise( vectorFunctionTools::filter( k, functions ), level, flag );
      }
   }

   void interpolate( const std::function< valueType( const hyteg::Point3D& ) >& expr, size_t level, DoFType flag = All ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); k++ )
      {
         compFunc_[k]->interpolate( expr, level, flag );
      }
   }

   void interpolate( valueType constant, size_t level, DoFType flag = All ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); k++ )
      {
         compFunc_[k]->interpolate( constant, level, flag );
      }
   }

   void interpolate( const std::vector< std::function< valueType( const hyteg::Point3D& ) > >& expr,
                     size_t                                                                    level,
                     DoFType                                                                   flag = All ) const
   {
      WALBERLA_ASSERT_GREATER_EQUAL( expr.size(), compFunc_.size() );
      WALBERLA_DEBUG_SECTION()
      {
         if ( expr.size() > compFunc_.size() )
         {
            WALBERLA_LOG_WARNING_ON_ROOT( "CSFVectorFunction::interpolate(): Ignoring excess std::functions!" );
         }
      }
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->interpolate( expr[k], level, flag );
      }
   }

   void interpolate( const std::vector< std::function< valueType( const Point3D& ) > >& expr,
                     uint_t                                                             level,
                     BoundaryUID                                                        boundaryUID ) const
   {
      WALBERLA_ASSERT_GREATER_EQUAL( expr.size(), compFunc_.size() );
      WALBERLA_DEBUG_SECTION()
      {
         if ( expr.size() > compFunc_.size() )
         {
            WALBERLA_LOG_WARNING_ON_ROOT( "CSFVectorFunction::interpolate(): Ignoring excess std::functions!" );
         }
      }
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->interpolate( expr[k], level, boundaryUID );
      }
   }

   void add( valueType scalar, size_t level, DoFType flag = All ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->add( scalar, level, flag );
      }
   }

   void add( const std::vector< valueType >                                           scalars,
             const std::vector< std::reference_wrapper< const VectorFunctionType > >& functions,
             size_t                                                                   level,
             DoFType                                                                  flag = All ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->add( scalars, vectorFunctionTools::filter( k, functions ), level, flag );
      }
   }

   void enableTiming( const std::shared_ptr< walberla::WcTimingTree >& timingTree )
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->enableTiming( timingTree );
      }
   }

   /// @name Bondary Conditions
   /// Methods for handling boundary conditions
   /// @{

   /// Returns a BoundaryCondition object, which is identical for all component functions
   BoundaryCondition getBoundaryCondition() const { return compFunc_[0]->getBoundaryCondition(); }

   /// Set boundary conditions, will be identical for all component functions
   void setBoundaryCondition( BoundaryCondition bc )
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->setBoundaryCondition( bc );
      }
   }

   template < typename OtherType >
   void copyBoundaryConditionFromFunction( const CSFVectorFunction< OtherType >& other )
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->setBoundaryCondition( vectorFunctionTools::filter( k, other ).getBoundaryCondition() );
      }
   }
   /// @}

   void swap( const VectorFunctionType& other, const uint_t& level, const DoFType& flag = All ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->swap( other[k], level, flag );
      }
   }

   void assign( const std::vector< valueType >&                                          scalars,
                const std::vector< std::reference_wrapper< const VectorFunctionType > >& functions,
                size_t                                                                   level,
                DoFType                                                                  flag = All ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->assign( scalars, vectorFunctionTools::filter( k, functions ), level, flag );
      }
   }

   valueType dotLocal( const VectorFunctionType& rhs, const uint_t level, const DoFType flag = All ) const
   {
      valueType sum = valueType( 0 );
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         sum += compFunc_[k]->dotLocal( rhs[k], level, flag );
      }
      return sum;
   }

   valueType dotGlobal( const VectorFunctionType& rhs, const uint_t level, const DoFType flag = All ) const
   {
      auto sum = dotLocal( rhs, level, flag );
      walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
      return sum;
   }

   valueType getMaxComponentMagnitude( uint_t level, DoFType flag, bool mpiReduce = true ) const
   {
      std::vector< valueType > values;

      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         values.push_back( compFunc_[k]->getMaxMagnitude( level, flag, false ) );
      }

      valueType localMax = *std::max_element( values.begin(), values.end() );

      valueType globalMax = localMax;
      if ( mpiReduce )
      {
         globalMax = walberla::mpi::allReduce( localMax, walberla::mpi::MAX );
      }

      return globalMax;
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
   void copyFrom( const VectorFunctionType&                      other,
                  const uint_t&                                  level,
                  const std::map< PrimitiveID::IDType, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID::IDType, uint_t >& otherPrimitiveIDsToRank ) const
   {
      for ( uint_t k = 0; k < compFunc_.size(); ++k )
      {
         compFunc_[k]->copyFrom( other[k], level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
      }
   }

 protected:
   const std::string                                     functionName_;
   std::vector< std::shared_ptr< VectorComponentType > > compFunc_;
};

} // namespace hyteg
