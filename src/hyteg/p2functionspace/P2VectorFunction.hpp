/*
 * Copyright (c) 2017-2020 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl, Andreas Wagner.
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

#include "hyteg/FunctionTraits.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"

namespace hyteg {

template < typename ValueType >
class P2VectorFunction
{
 public:
   using valueType = ValueType;

   using FunctionType = P2VectorFunction< ValueType >;

   using VectorComponentType = P2Function< ValueType >;

   using Tag = typename FunctionTrait< P2VectorFunction< ValueType > >::Tag;

   P2VectorFunction( const std::string&                         _name,
                     const std::shared_ptr< PrimitiveStorage >& storage,
                     size_t                                     minLevel,
                     size_t                                     maxLevel )
   : u( _name + "_u", storage, minLevel, maxLevel )
   , v( _name + "_v", storage, minLevel, maxLevel )
   , w( storage->hasGlobalCells() ? P2Function< ValueType >( _name + "_w", storage, minLevel, maxLevel ) :
                                    P2Function< ValueType >( _name + "_w_dummy", storage ) )

   , functionName_( _name )
   {}

   P2VectorFunction( const std::string&                         _name,
                     const std::shared_ptr< PrimitiveStorage >& storage,
                     size_t                                     minLevel,
                     size_t                                     maxLevel,
                     BoundaryCondition                          velocityBC )
   : u( _name + "_u", storage, minLevel, maxLevel, velocityBC )
   , v( _name + "_v", storage, minLevel, maxLevel, velocityBC )
   , w( storage->hasGlobalCells() ? P2Function< ValueType >( _name + "_w", storage, minLevel, maxLevel, velocityBC ) :
                                    P2Function< ValueType >( _name + "_w_dummy", storage ) )
   , functionName_( _name )
   {}

   std::shared_ptr< PrimitiveStorage > getStorage() const { return u.getStorage(); }

   bool isDummy() const { return false; }

   void interpolate( const std::function< ValueType( const hyteg::Point3D& ) >& expr, size_t level, DoFType flag = All ) const
   {
      u.interpolate( expr, level, flag );
      v.interpolate( expr, level, flag );
      w.interpolate( expr, level, flag );
   }

   void interpolate( const std::vector< std::function< ValueType( const hyteg::Point3D& ) > >& expr,
                     size_t                                                                    level,
                     DoFType                                                                   flag = All ) const
   {
      WALBERLA_ASSERT_GREATER( expr.size(), 0 );
      WALBERLA_ASSERT_LESS_EQUAL( expr.size(), 3 );

      for ( uint_t idx = 0; idx < expr.size(); ++idx )
      {
         component( idx ).interpolate( expr[idx], level, flag );
      }
   }

   void interpolate( const ValueType& constant, size_t level, DoFType flag = All ) const
   {
      u.interpolate( constant, level, flag );
      v.interpolate( constant, level, flag );
      w.interpolate( constant, level, flag );
   }

   void swap( const P2VectorFunction< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const
   {
      u.swap( other.u, level, flag );
      v.swap( other.v, level, flag );
      w.swap( other.w, level, flag );
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
   void copyFrom( const P2VectorFunction< ValueType >&           other,
                  const uint_t&                                  level,
                  const std::map< PrimitiveID::IDType, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID::IDType, uint_t >& otherPrimitiveIDsToRank ) const
   {
      u.copyFrom( other.u, level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
      v.copyFrom( other.v, level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
      w.copyFrom( other.w, level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
   }

   void assign( const std::vector< walberla::real_t >                                               scalars,
                const std::vector< std::reference_wrapper< const P2VectorFunction< ValueType > > >& functions,
                size_t                                                                              level,
                DoFType                                                                             flag = All ) const
   {
      std::vector< std::reference_wrapper< const P2Function< ValueType > > > functions_u;
      std::vector< std::reference_wrapper< const P2Function< ValueType > > > functions_v;
      std::vector< std::reference_wrapper< const P2Function< ValueType > > > functions_w;

      for ( const P2VectorFunction< ValueType >& function : functions )
      {
         functions_u.push_back( function.u );
         functions_v.push_back( function.v );
         functions_w.push_back( function.w );
      }

      u.assign( scalars, functions_u, level, flag );
      v.assign( scalars, functions_v, level, flag );
      w.assign( scalars, functions_w, level, flag );
   }

   void add( const std::vector< walberla::real_t >                                               scalars,
             const std::vector< std::reference_wrapper< const P2VectorFunction< ValueType > > >& functions,
             size_t                                                                              level,
             DoFType                                                                             flag = All ) const
   {
      std::vector< std::reference_wrapper< const P2Function< ValueType > > > functions_u;
      std::vector< std::reference_wrapper< const P2Function< ValueType > > > functions_v;
      std::vector< std::reference_wrapper< const P2Function< ValueType > > > functions_w;

      for ( const P2VectorFunction< ValueType >& function : functions )
      {
         functions_u.push_back( function.u );
         functions_v.push_back( function.v );
         functions_w.push_back( function.w );
      }

      u.add( scalars, functions_u, level, flag );
      v.add( scalars, functions_v, level, flag );
      w.add( scalars, functions_w, level, flag );
   }

   walberla::real_t dotLocal( const P2VectorFunction< ValueType >& rhs, const size_t level, const DoFType flag = All ) const
   {
      walberla::real_t sum = u.dotLocal( rhs.u, level, flag );
      sum += v.dotLocal( rhs.v, level, flag );
      sum += w.dotLocal( rhs.w, level, flag );
      return sum;
   }

   walberla::real_t dotGlobal( const P2VectorFunction< ValueType >& rhs, const size_t level, const DoFType flag = All ) const
   {
      auto sum = dotLocal( rhs, level, flag );
      walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
      return sum;
   }

   void enableTiming( const std::shared_ptr< walberla::WcTimingTree >& timingTree ) const
   {
      u.enableTiming( timingTree );
      v.enableTiming( timingTree );
      w.enableTiming( timingTree );
   }

   P2Function< ValueType >& component( uint_t idx )
   {
      switch ( idx )
      {
      case 0:
         return u;
      case 1:
         return v;
      case 2:
         return w;
      default:
         WALBERLA_ABORT( "Index out of range! Must be in {0,1,2}" );
      }
   }

   const P2Function< ValueType >& component( uint_t idx ) const
   {
      switch ( idx )
      {
      case 0:
         return u;
      case 1:
         return v;
      case 2:
         return w;
      default:
         WALBERLA_ABORT( "Index out of range! Must be in {0,1,2}" );
      }
   }

   const P2Function< ValueType >& operator[]( uint_t idx ) const
   {
      switch ( idx )
      {
      case 0:
         return u;
      case 1:
         return v;
      case 2:
         return w;
      default:
         WALBERLA_ABORT( "Index out of range! Must be in {0,1,2}" );
      }
   }

   P2Function< ValueType >& operator[]( uint_t idx )
   {
      switch ( idx )
      {
      case 0:
         return u;
      case 1:
         return v;
      case 2:
         return w;
      default:
         WALBERLA_ABORT( "Index out of range! Must be in {0,1,2}" );
      }
   }

   P2Function< ValueType > u;
   P2Function< ValueType > v;
   P2Function< ValueType > w;

   const std::string& getFunctionName() const { return functionName_; }

   uint_t getDimension() const { return component( 0 ).getStorage()->hasGlobalCells() ? 3 : 2; }

 private:
   const std::string functionName_;
};

} // namespace hyteg
