/*
 * Copyright (c) 2020-2023 Marcus Mohr.
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

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"
#include "hyteg/types/types.hpp"

namespace hyteg {

using walberla::uint_t;

// forward declaration of child for use in mother class
template < typename func_t >
class FunctionWrapper;

template < typename value_t >
class GenericFunction
{
 public:
   virtual ~GenericFunction(){};

   typedef value_t valueType;

   template < typename func_t >
   func_t& unwrap()
   {
      auto realMe = static_cast< FunctionWrapper< func_t >* >( this );
      return realMe->unwrap();
   }

   template < typename func_t >
   const func_t& unwrap() const
   {
      auto realMe = static_cast< const FunctionWrapper< func_t >* >( this );
      return realMe->unwrap();
   }

   /// Returns the dimension of the field represented by the function
   ///
   /// A function in HyTeG
   ///
   ///  Value  | Explanation
   /// :-----: | :-----------------
   ///     1   | a scalar field represented by a Function
   ///   d>1   | a vector field with components of dimension d represented by a VectorFunction
   ///     0   | a composite function such as e.g. from a Taylor-Hood discretisation, represented by a BlockFunction
   virtual uint_t getDimension() const = 0;

   virtual uint_t getMinLevel() const = 0;

   virtual uint_t getMaxLevel() const = 0;

   virtual const std::string& getFunctionName() const = 0;

   virtual functionTraits::FunctionKind getFunctionKind() const = 0;

   virtual uint_t getNumberOfLocalDoFs( uint_t level ) const = 0;

   virtual uint_t getNumberOfGlobalDoFs( uint_t          level,
                                         const MPI_Comm& communicator = walberla::mpi::MPIManager::instance()->comm(),
                                         const bool&     onRootOnly   = false ) const = 0;

   virtual std::shared_ptr< PrimitiveStorage > getStorage() const = 0;

   virtual void multElementwise( const std::vector< std::reference_wrapper< const GenericFunction< value_t > > >& functions,
                                 uint_t                                                                           level,
                                 DoFType flag = All ) const = 0;

   virtual void interpolate( value_t constant, uint_t level, DoFType flag = All ) const = 0;

   virtual value_t
       dotGlobal( const GenericFunction< value_t >& secondOp, const uint_t level, const DoFType flag = All ) const = 0;

   virtual value_t dotLocal( const GenericFunction< value_t >& secondOp, uint_t level, DoFType flag = All ) const = 0;

   virtual void enableTiming( const std::shared_ptr< walberla::WcTimingTree >& timingTree ) = 0;

   virtual void setBoundaryCondition( BoundaryCondition bc ) = 0;

   virtual BoundaryCondition getBoundaryCondition() const = 0;

   virtual void add( const value_t scalar, uint_t level, DoFType flag = All ) const = 0;

   virtual void add( const std::vector< value_t >                                                     scalars,
                     const std::vector< std::reference_wrapper< const GenericFunction< value_t > > >& functions,
                     uint_t                                                                           level,
                     DoFType                                                                          flag = All ) const = 0;

   virtual void assign( const std::vector< value_t >                                                     scalars,
                        const std::vector< std::reference_wrapper< const GenericFunction< value_t > > >& functions,
                        uint_t                                                                           level,
                        DoFType                                                                          flag = All ) const = 0;

   virtual void
       interpolate( const std::function< value_t( const hyteg::Point3D& ) >& expr, uint_t level, DoFType flag = All ) const = 0;

   virtual void interpolate( const std::vector< std::function< value_t( const hyteg::Point3D& ) > >& expr,
                             uint_t                                                                  level,
                             DoFType                                                                 flag = All ) const = 0;

   virtual void swap( const GenericFunction< value_t >& other, const uint_t& level, const DoFType& flag = All ) const = 0;

   virtual void copyFrom( const GenericFunction< value_t >&      other,
                          const uint_t&                          level,
                          const std::map< PrimitiveID, uint_t >& localPrimitiveIDsToRank,
                          const std::map< PrimitiveID, uint_t >& otherPrimitiveIDsToRank ) const = 0;

   // currently there seems no way around providing only a few type alternatives
   // fixing this would require a restructuring of GenericFunction and FunctionWrapper as there are no templated virtual functions
   virtual void copyBoundaryConditionFromFunction( const GenericFunction< double >& other ) = 0;
   virtual void copyBoundaryConditionFromFunction( const GenericFunction< float >& other ) = 0;
   virtual void copyBoundaryConditionFromFunction( const GenericFunction< idx_t >& other ) = 0;

   virtual void enumerate( uint_t level ) const = 0;

   virtual void enumerate( uint_t level, value_t& offset ) const = 0;

   /// conversion to/from linear algebra representation
   /// @{
   virtual void toVector( const GenericFunction< idx_t >&       numerator,
                          const std::shared_ptr< VectorProxy >& vec,
                          uint_t                                level,
                          DoFType                               flag ) const = 0;

   virtual void fromVector( const GenericFunction< idx_t >&       numerator,
                            const std::shared_ptr< VectorProxy >& vec,
                            uint_t                                level,
                            DoFType                               flag ) const = 0;
   /// @}
};

} // namespace hyteg
