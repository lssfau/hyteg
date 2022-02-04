/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"

namespace hyteg {

template < typename ValueType >
class P2P2StokesFunction
{
 public:
   using valueType = ValueType;

   template < typename VType >
   using FunctionType = P2P2StokesFunction< VType >;

   using VelocityFunction_T = P2VectorFunction< ValueType >;
   using PressureFunction_T = P2Function< ValueType >;

   using Tag = typename FunctionTrait< P2P2StokesFunction< ValueType > >::Tag;

   P2P2StokesFunction( const std::string&                         _name,
                       const std::shared_ptr< PrimitiveStorage >& storage,
                       size_t                                     minLevel,
                       size_t                                     maxLevel )
   : uvw( _name + "_vector", storage, minLevel, maxLevel )
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

   void swap( const P2P2StokesFunction< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const
   {
      uvw.swap( other.uvw, level, flag );
      p.swap( other.p, level, flag );
   }

   void assign( const std::vector< walberla::real_t >                                                 scalars,
                const std::vector< std::reference_wrapper< const P2P2StokesFunction< ValueType > > >& functions,
                size_t                                                                                level,
                DoFType                                                                               flag = All ) const
   {
      std::vector< std::reference_wrapper< const P2VectorFunction< ValueType > > > functions_uvw;
      std::vector< std::reference_wrapper< const P2Function< ValueType > > >       functions_p;

      for ( const P2P2StokesFunction< ValueType >& function : functions )
      {
         functions_uvw.push_back( function.uvw );
         functions_p.push_back( function.p );
      }

      uvw.assign( scalars, functions_uvw, level, flag );
      p.assign( scalars, functions_p, level, flag );
   }

   void add( const std::vector< walberla::real_t >                                                 scalars,
             const std::vector< std::reference_wrapper< const P2P2StokesFunction< ValueType > > >& functions,
             size_t                                                                                level,
             DoFType                                                                               flag = All ) const
   {
      std::vector< std::reference_wrapper< const P2VectorFunction< ValueType > > > functions_uvw;
      std::vector< std::reference_wrapper< const P2Function< ValueType > > >       functions_p;

      for ( const P2P2StokesFunction< ValueType >& function : functions )
      {
         functions_uvw.push_back( function.uvw );
         functions_p.push_back( function.p );
      }

      uvw.add( scalars, functions_uvw, level, flag );
      p.add( scalars, functions_p, level, flag );
   }

   walberla::real_t dotGlobal( const P2P2StokesFunction< ValueType >& rhs, const uint_t level, const DoFType flag = All ) const
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

   /// conversion to/from linear algebra representation
   /// @{
   void toVector( const P2P2StokesFunction< idx_t >&      numerator,
                  const std::shared_ptr< VectorProxy >& vec,
                  uint_t                                level,
                  DoFType                               flag ) const
   {
      for ( uint_t k = 0; k < uvw.getDimension(); ++k )
      {
         uvw[k].toVector( numerator.uvw[k], vec, level, flag );
      }
      p.toVector( numerator.p, vec, level, flag );
   }

   void fromVector( const P2P2StokesFunction< idx_t >&      numerator,
                    const std::shared_ptr< VectorProxy >& vec,
                    uint_t                                level,
                    DoFType                               flag ) const
   {
      for ( uint_t k = 0; k < uvw.getDimension(); ++k )
      {
         uvw[k].fromVector( numerator.uvw[k], vec, level, flag );
      }
      p.fromVector( numerator.p, vec, level, flag );
   };
   /// @}

   P2VectorFunction< ValueType > uvw;
   P2Function< ValueType >       p;
};

} // namespace hyteg
