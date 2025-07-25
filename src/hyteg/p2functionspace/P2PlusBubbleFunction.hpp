/*
 * Copyright (c) 2017-2024 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/volumedofspace/VolumeDoFFunction.hpp"

namespace hyteg {

using walberla::real_c;

/// Class to represent a continuous Lagrange function  of second order enhanced by a bubble
template < typename ValueType >
class P2PlusBubbleFunction final : public Function< P2PlusBubbleFunction< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = P2PlusBubbleFunction< VType >;

   P2PlusBubbleFunction( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage );

   P2PlusBubbleFunction( const std::string&                         name,
                         const std::shared_ptr< PrimitiveStorage >& storage,
                         uint_t                                     minLevel,
                         uint_t                                     maxLevel );

   P2PlusBubbleFunction( const std::string&                         name,
                         const std::shared_ptr< PrimitiveStorage >& storage,
                         uint_t                                     minLevel,
                         uint_t                                     maxLevel,
                         BoundaryCondition                          boundaryCondition );

   uint_t getDimension() const override final { return 1; }

   P2PlusBubbleFunction< ValueType >& operator[]( uint_t idx )
   {
      WALBERLA_ASSERT( idx == 0 );
      WALBERLA_UNUSED( idx );
      return *this;
   }

   const P2PlusBubbleFunction< ValueType >& operator[]( uint_t idx ) const
   {
      WALBERLA_ASSERT( idx == 0 );
      WALBERLA_UNUSED( idx );
      return *this;
   }

   inline const vertexdof::VertexDoFFunction< ValueType >&      getVertexDoFFunction() const { return vertexDoFFunction_; }
   inline const EdgeDoFFunction< ValueType >&                   getEdgeDoFFunction() const { return edgeDoFFunction_; }
   inline const volumedofspace::VolumeDoFFunction< ValueType >& getVolumeDoFFunction() const { return volumeDoFFunction_; }

   template < typename SenderType, typename ReceiverType >
   void startCommunication( const uint_t& level ) const
   {
      WALBERLA_ABORT( "Member function not bubble ready!" );
      vertexDoFFunction_.template startCommunication< SenderType, ReceiverType >( level );
      edgeDoFFunction_.template startCommunication< SenderType, ReceiverType >( level );
   }

   template < typename SenderType, typename ReceiverType >
   void endCommunication( const uint_t& level ) const
   {
      WALBERLA_ABORT( "Member function not bubble ready!" );
      vertexDoFFunction_.template endCommunication< SenderType, ReceiverType >( level );
      edgeDoFFunction_.template endCommunication< SenderType, ReceiverType >( level );
   }

   template < typename SenderType, typename ReceiverType >
   void communicate( const uint_t& level ) const
   {
      vertexDoFFunction_.template communicate< SenderType, ReceiverType >( level );
      edgeDoFFunction_.template communicate< SenderType, ReceiverType >( level );
   }

   /// \brief Evaluate finite element function at a specific coordinates.
   ///
   /// In a parallel setting, the specified coordinate might not lie in the local subdomain.
   ///
   /// Evaluation is performed in two steps:
   ///
   ///   1. For all volume primitives of the local subdomain:
   ///      If a point-tet (point-triangle in 2D) inclusion test succeeds, the function returns true
   ///      and the finite-element function is evaluated.
   ///
   ///   2. Skipped, if radius is negative.
   ///      For all volume primitives of the local subdomain:
   ///      A sphere-tet (circle-triangle) intersection
   ///      test is performed, if successful returns true, and the finite-element function is extrapolated
   ///      to the specified coordinate and evaluated, depending on the radius, this might introduce (large) errors.
   ///
   /// If both tests fail, this function returns false, and no evaluation is performed (i.e. the returned, evaluated
   /// value is not set to anything meaningful).
   ///
   /// Note that two parallel processes that return true, may return _different_ values.
   ///
   /// No communication is performed in this function.
   /// -> Does not need to be called collectively.
   /// -> Different values are returned on each process.
   ///
   /// \param physicalCoords coordinates in physical domain where the function is to be evaluated
   /// \param level refinement level
   /// \param value function value at the coordinate if search was successful
   /// \param searchToleranceRadius radius of the sphere (circle) for the second search phase, skipped if negative
   /// \return true if the function was evaluated successfully, false otherwise
   ///
   bool evaluate( const Point3D& physicalCoords,
                  uint_t         level,
                  ValueType&     value,
                  real_t         searchToleranceRadius = real_c( 1e-05 ) ) const;

   inline void evaluateGradient( const Point3D& physicalCoords, uint_t level, Point3D& gradient ) const;

   /// @name Member functions for interpolation using BoundaryUID flags
   ///@{
   void interpolate( ValueType constant, uint_t level, BoundaryUID boundaryUID ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, BoundaryUID boundaryUID ) const;
   ///@}

   /// @name Member functions for interpolation using DoFType flags
   ///@{
   void interpolate( ValueType constant, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::vector< std::function< ValueType( const Point3D& ) > >& expr,
                     uint_t                                                             level,
                     DoFType                                                            flag = All ) const
   {
      WALBERLA_ABORT( "Member function not bubble ready!" );
      WALBERLA_ASSERT_EQUAL( expr.size(), 1 );
      this->interpolate( expr[0], level, flag );
   };

   void interpolate( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >&    expr,
                     const std::vector< std::reference_wrapper< const P2PlusBubbleFunction< ValueType > > >& srcFunctions,
                     uint_t                                                                                  level,
                     DoFType                                                                                 flag = All ) const;
   ///@}

   /// Set all function DoFs to zero including the ones in the halos
   void setToZero( const uint_t level ) const override final;

   void swap( const P2PlusBubbleFunction< ValueType >& other, const uint_t& level, const DoFType& dofType = All ) const;

   /// \brief Copies all values function data from other to this.
   ///
   /// This method can be used safely if the other function is located on a different PrimitiveStorage.
   void copyFrom( const P2PlusBubbleFunction< ValueType >& other, const uint_t& level ) const;

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
   void copyFrom( const P2PlusBubbleFunction< ValueType >& other,
                  const uint_t&                            level,
                  const std::map< PrimitiveID, uint_t >&   localPrimitiveIDsToRank,
                  const std::map< PrimitiveID, uint_t >&   otherPrimitiveIDsToRank ) const;

   void assign( const std::vector< ValueType >&                                                         scalars,
                const std::vector< std::reference_wrapper< const P2PlusBubbleFunction< ValueType > > >& functions,
                uint_t                                                                                  level,
                DoFType                                                                                 flag = All ) const;

   void add( ValueType scalar, uint_t level, DoFType flag = All ) const;

   void add( const std::vector< ValueType >&                                                         scalars,
             const std::vector< std::reference_wrapper< const P2PlusBubbleFunction< ValueType > > >& functions,
             uint_t                                                                                  level,
             DoFType                                                                                 flag = All ) const;

   /// Compute the product of several functions in an elementwise fashion
   ///
   /// The method takes as input a collection of functions. These are multiplied together in an elementwise fashion.
   /// The latter is to be understood not in a FE context, but in the sense of element-wise operators in matrix/array
   /// oriented languages, i.e. the product is a function of the same type as the inputs and its DoFs are formed as
   /// product of the corresponding DoFs of the input functions. The result is stored in the function object on which
   /// the method is invoked, overwritting its contents. It is safe, if the destination function is part of the product.
   ///
   /// \param functions  the functions forming the product
   /// \param level      level on which the multiplication should be computed
   /// \param flag       marks those primitives which are partaking in the computation of the product
   void multElementwise( const std::vector< std::reference_wrapper< const P2PlusBubbleFunction< ValueType > > >& functions,
                         uint_t                                                                                  level,
                         DoFType flag = All ) const;

   /// Replace values of the function by their inverses in an elementwise fashion
   void invertElementwise( uint_t level, DoFType flag = All, bool workOnHalos = false ) const;

   ValueType dotGlobal( const P2PlusBubbleFunction< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const;

   ValueType dotLocal( const P2PlusBubbleFunction< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const;

   ValueType sumGlobal( uint_t level, const DoFType& flag = All, const bool& absolute = false ) const;

   ValueType sumLocal( uint_t level, const DoFType& flag = All, const bool& absolute = false ) const;

   void
       prolongateP1ToP2( const hyteg::P1Function< ValueType >& p1Function, const uint_t& level, const DoFType& flag = All ) const;

   void restrictP2ToP1( const P1Function< ValueType >& p1Function, const uint_t& level, const DoFType& flag = All ) const;

   void restrictInjection( uint_t sourceLevel, DoFType flag = All ) const;

   ValueType getMaxDoFValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;

   ValueType getMaxDoFMagnitude( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;

   ValueType getMinDoFValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;

   BoundaryCondition getBoundaryCondition() const;
   void              setBoundaryCondition( BoundaryCondition bc );

   template < typename OtherFunctionValueType >
   void copyBoundaryConditionFromFunction( const P2PlusBubbleFunction< OtherFunctionValueType >& other )
   {
      setBoundaryCondition( other.getBoundaryCondition() );
   }

   void enumerate( uint_t level ) const;

   void enumerate( uint_t level, ValueType& offset ) const;

   void setLocalCommunicationMode( const communication::BufferedCommunicator::LocalCommunicationMode& localCommMode );

   /// @name Member functions for conversion to/from linear algebra representation
   ///@{
   void toVector( const P2PlusBubbleFunction< idx_t >&  numerator,
                  const std::shared_ptr< VectorProxy >& vec,
                  uint_t                                level,
                  DoFType                               flag ) const
   {
      this->getVertexDoFFunction().toVector( numerator.getVertexDoFFunction(), vec, level, flag );
      this->getEdgeDoFFunction().toVector( numerator.getEdgeDoFFunction(), vec, level, flag );
      this->getVolumeDoFFunction().toVector( numerator.getVolumeDoFFunction(), vec, level, flag );
   }

   void fromVector( const P2PlusBubbleFunction< idx_t >&  numerator,
                    const std::shared_ptr< VectorProxy >& vec,
                    uint_t                                level,
                    DoFType                               flag ) const
   {
      this->getVertexDoFFunction().fromVector( numerator.getVertexDoFFunction(), vec, level, flag );
      this->getEdgeDoFFunction().fromVector( numerator.getEdgeDoFFunction(), vec, level, flag );
      this->getVolumeDoFFunction().fromVector( numerator.getVolumeDoFFunction(), vec, level, flag );
   };
   ///@}

   /// \brief Returns the number of DoFs. Performs global reduction, must be called collectively.
   ///
   /// \param level        refinement level
   /// \param communicator if required, a custom communicator can be passed
   /// \param onRootOnly   if true, the result is only returned on the root process
   /// \return
   uint_t getNumberOfGlobalDoFs( uint_t          level,
                                 const MPI_Comm& communicator = walberla::mpi::MPIManager::instance()->comm(),
                                 const bool&     onRootOnly   = false ) const
   {
      return numberOfGlobalDoFs< P2PlusBubbleFunctionTag >( *this->storage_, level, communicator, onRootOnly );
   }

 private:
   using Function< P2PlusBubbleFunction< ValueType > >::communicators_;

   vertexdof::VertexDoFFunction< ValueType >      vertexDoFFunction_;
   EdgeDoFFunction< ValueType >                   edgeDoFFunction_;
   volumedofspace::VolumeDoFFunction< ValueType > volumeDoFFunction_;
};

extern template class P2PlusBubbleFunction< double >;
extern template class P2PlusBubbleFunction< float >;
extern template class P2PlusBubbleFunction< int32_t >;
extern template class P2PlusBubbleFunction< idx_t >;

} //namespace hyteg
