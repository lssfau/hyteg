/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

template < typename ValueType >
class P1DGEFunction final : public Function< P1DGEFunction< ValueType > >
{
 public:
   P1DGEFunction( const std::string&                         name,
                  const std::shared_ptr< PrimitiveStorage >& storage,
                  uint_t                                     minLevel,
                  uint_t                                     maxLevel,
                  BoundaryCondition                          boundaryCondition = BoundaryCondition::create0123BC() );

   std::shared_ptr< P1VectorFunction< ValueType > > getConformingPart() const { return u_conforming_; }

   std::shared_ptr< dg::DGFunction< ValueType > > getDiscontinuousPart() const { return u_discontinuous_; }

   void setBoundaryCondition( BoundaryCondition bc )
   {
      u_conforming_->setBoundaryCondition( bc );
      u_discontinuous_->setBoundaryCondition( bc );
   }

   [[nodiscard]] BoundaryCondition getBoundaryCondition() const { return u_conforming_->getBoundaryCondition(); }

   template < typename SenderType, typename ReceiverType >
   void communicate( const uint_t& level ) const
   {
      u_conforming_->communicate( level );
      u_discontinuous_->communicate( level );
   }

   void add( const ValueType scalar, uint_t level, DoFType flag = All ) const
   {
      u_conforming_->add( scalar, level, flag );
      u_discontinuous_->add( scalar, level, flag );
   };

   void add( const std::vector< ValueType >                                                   scalars,
             const std::vector< std::reference_wrapper< const P1DGEFunction< ValueType > > >& functions,
             uint_t                                                                           level,
             DoFType                                                                          flag = All ) const
   {
      u_conforming_->add( scalars, functions, level, flag );
      u_discontinuous_->add( scalars, functions, level, flag );
   };

   void multElementwise( const std::vector< std::reference_wrapper< const P1DGEFunction< ValueType > > >& functions,
                         uint_t                                                                           level,
                         DoFType                                                                          flag = All ) const
   {
      u_conforming_->multElementwise( functions, level, flag );
      u_discontinuous_->multElementwise( functions, level, flag );
   }

   void interpolate( ValueType constant, uint_t level, DoFType dofType = All ) const
   {
      WALBERLA_UNUSED( dofType );
      WALBERLA_ABORT( "Not implemented" );
   }

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType dofType = All ) const
   {
      WALBERLA_UNUSED( dofType );
      WALBERLA_ABORT( "Not implemented" );
   }

   void interpolate( const std::vector< std::function< ValueType( const hyteg::Point3D& ) > >& expressions,
                     uint_t                                                                    level,
                     DoFType                                                                   flag = All ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   };

   void swap( const P1DGEFunction< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const
   {
      u_conforming_->swap( other.u_conforming_, level, flag );
      u_discontinuous_->swap( other.u_discontinuous_, level, flag );
   };

   void copyFrom( const P0Function< ValueType >&                 other,
                  const uint_t&                                  level,
                  const std::map< PrimitiveID::IDType, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID::IDType, uint_t >& otherPrimitiveIDsToRank ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   };

   bool evaluate( const Point3D& coordinates, uint_t level, ValueType& value, real_t searchToleranceRadius = 1e-05 ) const
   {
      WALBERLA_ABORT( "Not implemented." );
   }

   void assign( const std::vector< ValueType >&                                                  scalars,
                const std::vector< std::reference_wrapper< const P1DGEFunction< ValueType > > >& functions,
                uint_t                                                                           level,
                DoFType                                                                          flag = All ) const
   {
      std::vector< std::reference_wrapper< const dg::DGFunction< ValueType > > >   dg_list;
      std::vector< std::reference_wrapper< const P1VectorFunction< ValueType > > > p1_list;

      for ( auto f : functions )
      {
         dg_list.push_back( *f.getDiscontinuousPart() );
         p1_list.push_back( *f.getConformingPart() );
      }
      u_conforming_->assign( scalars, dg_list, level, flag );
      u_discontinuous_->assign( scalars, p1_list, level, flag );
   }

   ValueType dotGlobal( const P1DGEFunction< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const
   {
      real_t result = 0;
      result += u_conforming_->dotGlobal( *rhs.getConformingPart(), level, flag );
      result += u_discontinuous_->dotGlobal( *rhs.getDiscontinuousPart(), level, flag );
      return result;
   }

   ValueType dotLocal( const P1DGEFunction< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const
   {
      real_t result = 0;
      result += u_conforming_->dotLlobal( *rhs.getConformingPart(), level, flag );
      result += u_discontinuous_->dotLlobal( *rhs.getDiscontinuousPart(), level, flag );
      return result;
   }

   void enumerate( uint_t level ) const { enumerate( level, 0 ); }

   void enumerate( uint_t level, ValueType& offset ) const
   {
      u_conforming_->enumerate( level );
      offset += u_conforming_->getNumberOfGlobalDoFs( level );
      u_discontinuous_->enumerate( level, offset );
   }

   /// conversion to/from linear algebra representation
   /// @{
   void toVector( const P1DGEFunction< idx_t >&         numerator,
                  const std::shared_ptr< VectorProxy >& vec,
                  uint_t                                level,
                  DoFType                               flag ) const
   {
      u_conforming_->toVector( *numerator.getConformingPart(), vec, level, flag );
      u_discontinuous_->toVector( *numerator.getDiscontinuousPart(), vec, level, flag );
   }

   void fromVector( const P1DGEFunction< idx_t >&         numerator,
                    const std::shared_ptr< VectorProxy >& vec,
                    uint_t                                level,
                    DoFType                               flag ) const
   {
      u_conforming_->fromVector( *numerator.getConformingPart(), vec, level, flag );
      u_discontinuous_->fromVector( *numerator.getDiscontinuousPart(), vec, level, flag );
   }
   /// @}

   [[nodiscard]] uint_t getNumberOfLocalDoFs( uint_t level ) const
   {
      return u_conforming_->getNumberOfLocalDoFs( level ) + u_discontinuous_->getNumberOfLocalDoFs( level );
   }

   [[nodiscard]] uint_t getNumberOfGlobalDoFs( uint_t          level,
                                               const MPI_Comm& communicator = walberla::mpi::MPIManager::instance()->comm(),
                                               const bool&     onRootOnly   = false ) const
   {
      return u_conforming_->getNumberOfGlobalDoFs( level ) + u_discontinuous_->getNumberOfGlobalDoFs( level );
   }

 protected:
   std::shared_ptr< dg::DGBasisLinearLagrange_Example > basis_;

   std::shared_ptr< P1VectorFunction< ValueType > > u_conforming_;
   std::shared_ptr< dg::DGFunction< ValueType > >   u_discontinuous_;
};

template < typename ValueType >
P1DGEFunction< ValueType >::P1DGEFunction( const std::string&                         name,
                                           const std::shared_ptr< PrimitiveStorage >& storage,
                                           uint_t                                     minLevel,
                                           uint_t                                     maxLevel,
                                           BoundaryCondition                          boundaryCondition )
: Function< P1DGEFunction< ValueType > >( name, storage, minLevel, maxLevel )
, basis_{ std::make_shared< dg::DGBasisLinearLagrange_Example >() }
, u_conforming_{ std::make_shared< P1VectorFunction< ValueType > >( name, storage, minLevel, maxLevel, boundaryCondition ) }
, u_discontinuous_{
      std::make_shared< dg::DGFunction< ValueType > >( name, storage, minLevel, maxLevel, basis_, 0, boundaryCondition ) }
{}

} // namespace hyteg
