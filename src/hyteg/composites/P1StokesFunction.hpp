/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"

namespace hyteg {

template < typename ValueType >
class P1StokesFunction
{
 public:
   using valueType = ValueType;

   template < typename VType >
   using FunctionType = P1StokesFunction< VType >;

   using VelocityFunction_T = P1VectorFunction< ValueType >;
   using PressureFunction_T = P1Function< ValueType >;

   using Tag = typename FunctionTrait< P1StokesFunction< ValueType > >::Tag;

   P1StokesFunction( const std::string&                         _name,
                     const std::shared_ptr< PrimitiveStorage >& storage,
                     size_t                                     minLevel,
                     size_t                                     maxLevel )
   : P1StokesFunction( _name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
   {}

   P1StokesFunction( const std::string&                         _name,
                     const std::shared_ptr< PrimitiveStorage >& storage,
                     size_t                                     minLevel,
                     size_t                                     maxLevel,
                     BoundaryCondition                          velocityBC )
   : uvw( _name + "_vector", storage, minLevel, maxLevel, velocityBC )
   , p( _name + "_p", storage, minLevel, maxLevel, BoundaryCondition::createAllInnerBC() )
   {}

   std::shared_ptr< PrimitiveStorage > getStorage() const { return uvw.getStorage(); }

   bool isDummy() const { return false; }

   void interpolate( const std::function< real_t( const hyteg::Point3D& ) >& expr, size_t level, DoFType flag = All ) const
   {
      uvw.interpolate( expr, level, flag );
      p.interpolate( expr, level, flag );
   }

   void interpolate( const real_t& constant, size_t level, DoFType flag = All ) const
   {
      uvw.interpolate( constant, level, flag );
      p.interpolate( constant, level, flag );
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
   void copyFrom( const P1StokesFunction< ValueType >&           other,
                  const uint_t&                                  level,
                  const std::map< PrimitiveID::IDType, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID::IDType, uint_t >& otherPrimitiveIDsToRank ) const
   {
      uvw.copyFrom( other.uvw, level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
      p.copyFrom( other.p, level, localPrimitiveIDsToRank, otherPrimitiveIDsToRank );
   }

   void swap( const P1StokesFunction< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const
   {
      uvw.swap( other.uvw, level, flag );
      p.swap( other.p, level, flag );
   }

   void assign( const std::vector< walberla::real_t >                                               scalars,
                const std::vector< std::reference_wrapper< const P1StokesFunction< ValueType > > >& functions,
                size_t                                                                              level,
                DoFType                                                                             flag = All ) const
   {
      std::vector< std::reference_wrapper< const P1VectorFunction< ValueType > > > functions_uvw;
      std::vector< std::reference_wrapper< const P1Function< ValueType > > >       functions_p;

      for ( const P1StokesFunction< ValueType >& function : functions )
      {
         functions_uvw.push_back( function.uvw );
         functions_p.push_back( function.p );
      }

      uvw.assign( scalars, functions_uvw, level, flag );
      p.assign( scalars, functions_p, level, flag );
   }

   void add( real_t scalar, size_t level, DoFType flag = All ) const
   {
      uvw.add( scalar, level, flag );
      p.add( scalar, level, flag );
   }

   void add( const std::vector< walberla::real_t >                                               scalars,
             const std::vector< std::reference_wrapper< const P1StokesFunction< ValueType > > >& functions,
             size_t                                                                              level,
             DoFType                                                                             flag = All ) const
   {
      std::vector< std::reference_wrapper< const P1VectorFunction< ValueType > > > functions_uvw;
      std::vector< std::reference_wrapper< const P1Function< ValueType > > >       functions_p;

      for ( const P1StokesFunction< ValueType >& function : functions )
      {
         functions_uvw.push_back( function.uvw );
         functions_p.push_back( function.p );
      }

      uvw.add( scalars, functions_uvw, level, flag );
      p.add( scalars, functions_p, level, flag );
   }

   walberla::real_t dotGlobal( const P1StokesFunction< ValueType >& rhs, const uint_t level, const DoFType flag = All ) const
   {
      walberla::real_t sum = uvw.dotLocal( rhs.uvw, level, flag );
      sum += p.dotLocal( rhs.p, level, flag );
      walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
      return sum;
   }

   void enableTiming( const std::shared_ptr< walberla::WcTimingTree >& timingTree )
   {
      uvw.enableTiming( timingTree );
      p.enableTiming( timingTree );
   }

   void enumerate( uint_t level )
   {
      uint_t counterVertexDoFs = hyteg::numberOfLocalDoFs< Tag >( *( getStorage() ), level );

      std::vector< uint_t > vertexDoFsPerRank = walberla::mpi::allGather( counterVertexDoFs );

      ValueType offset = 0;

      for ( uint_t i = 0; i < uint_c( walberla::MPIManager::instance()->rank() ); ++i )
      {
         offset += static_cast< ValueType >( vertexDoFsPerRank[i] );
      }

      for ( uint_t k = 0; k < uvw.getDimension(); k++ )
      {
         uvw[k].enumerate( level, offset );
      }
      p.enumerate( level, offset );
   }

   BoundaryCondition getVelocityBoundaryCondition() const { return uvw.getBoundaryCondition(); }

   BoundaryCondition getPressureBoundaryCondition() const { return p.getBoundaryCondition(); }

   void setVelocityBoundaryCondition( BoundaryCondition bc ) { uvw.setBoundaryCondition( bc ); }

   void setPressureBoundaryCondition( BoundaryCondition bc ) { p.setBoundaryCondition( bc ); }

   template < typename OtherFunctionValueType >
   void copyBoundaryConditionFromFunction( const P1StokesFunction< OtherFunctionValueType >& other )
   {
      setVelocityBoundaryCondition( other.getVelocityBoundaryCondition() );
      setPressureBoundaryCondition( other.getPressureBoundaryCondition() );
   }

   P1VectorFunction< ValueType > uvw;
   P1Function< ValueType >       p;
};

inline unsigned long long p1p1localFunctionMemorySize( const uint_t& level, const std::shared_ptr< PrimitiveStorage >& storage )
{
   if ( storage->hasGlobalCells() )
   {
      return 4 * vertexDoFLocalFunctionMemorySize( level, storage );
   }
   else
   {
      return 3 * vertexDoFLocalFunctionMemorySize( level, storage );
   }
}

inline unsigned long long p1p1globalFunctionMemorySize( const uint_t& level, const std::shared_ptr< PrimitiveStorage >& storage )
{
   const auto memLocal  = p1p1localFunctionMemorySize( level, storage );
   const auto memGlobal = walberla::mpi::allReduce( memLocal, walberla::mpi::SUM );
   return memGlobal;
}

} // namespace hyteg
