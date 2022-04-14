/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "core/mpi/all.h"

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

namespace hyteg {

using walberla::real_c;

template < typename ValueType >
class FunctionMemory;
template < typename DataType, typename PrimitiveType >
class PrimitiveDataID;

class PrimitiveStorage;
class Vertex;
class Edge;
class Face;
class Cell;

namespace edgedof {
///@name Size Functions
///@{

uint_t edgeDoFMacroVertexFunctionMemorySize( const uint_t& level, const Primitive& primitive );

uint_t edgeDoFMacroEdgeFunctionMemorySize( const uint_t& level, const Primitive& primitive );

uint_t edgeDoFMacroFaceFunctionMemorySize( const uint_t& level, const Primitive& primitive );

uint_t edgeDoFMacroCellFunctionMemorySize( const uint_t& level, const Primitive& primitive );

unsigned long long edgeDoFLocalFunctionMemorySize( const uint_t & level, const std::shared_ptr< PrimitiveStorage > & storage );

unsigned long long edgeDoFGlobalFunctionMemorySize( const uint_t & level, const std::shared_ptr< PrimitiveStorage > & storage );

///@}

} // namespace edgedof

template < typename ValueType >
class EdgeDoFFunction final : public Function< EdgeDoFFunction< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = EdgeDoFFunction< VType >;

   EdgeDoFFunction( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage );

   EdgeDoFFunction( const std::string&                         name,
                    const std::shared_ptr< PrimitiveStorage >& storage,
                    const uint_t&                              minLevel,
                    const uint_t&                              maxLevel );

   EdgeDoFFunction( const std::string&                         name,
                    const std::shared_ptr< PrimitiveStorage >& storage,
                    const uint_t&                              minLevel,
                    const uint_t&                              maxLevel,
                    const BoundaryCondition&                   boundaryCondition );

   bool hasMemoryAllocated( const uint_t & level, const Vertex & vertex ) const;
   bool hasMemoryAllocated( const uint_t & level, const Edge & edge ) const;
   bool hasMemoryAllocated( const uint_t & level, const Face & face ) const;
   bool hasMemoryAllocated( const uint_t & level, const Cell & cell ) const;

   void allocateMemory( const uint_t & level, const Vertex & vertex );
   void allocateMemory( const uint_t & level, const Edge & edge );
   void allocateMemory( const uint_t & level, const Face & face );
   void allocateMemory( const uint_t & level, const Cell & cell );

   void deleteMemory( const uint_t & level, const Vertex & vertex );
   void deleteMemory( const uint_t & level, const Edge & edge );
   void deleteMemory( const uint_t & level, const Face & face );
   void deleteMemory( const uint_t & level, const Cell & cell );

   void swap( const EdgeDoFFunction< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const;

   /// \brief Copies all values function data from other to this.
   ///
   /// This method can be used safely if the other function is located on a different PrimitiveStorage.
   void copyFrom( const EdgeDoFFunction< ValueType >& other, const uint_t& level ) const;

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
   void copyFrom( const EdgeDoFFunction< ValueType >&            other,
                  const uint_t&                                  level,
                  const std::map< PrimitiveID::IDType, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID::IDType, uint_t >& otherPrimitiveIDsToRank ) const;

   void assign( const std::vector< ValueType >&                                                    scalars,
                const std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >& functions,
                uint_t                                                                             level,
                DoFType                                                                            flag = All ) const;

   void add( ValueType scalar, uint_t level, DoFType flag = All ) const;

   void add( const std::vector< ValueType >&                                                    scalars,
             const std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >& functions,
             uint_t                                                                             level,
             DoFType                                                                            flag = All ) const;

   /// @name Member functions for interpolation using BoundaryUID flags
   //@{
   void interpolate( ValueType constant, uint_t level, BoundaryUID boundaryUID ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, BoundaryUID boundaryUID ) const;

   void interpolateExtended( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                             const std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >&   srcFunctions,
                             uint_t                                                                               level,
                             BoundaryUID boundaryUID ) const;
   //@}

   /// @name Member functions for interpolation using DoFType flags
   //@{
   void interpolate( ValueType constant, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType flag = All ) const;

   void interpolateExtended( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                             const std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >&   srcFunctions,
                             uint_t                                                                               level,
                             DoFType flag = All ) const;

   void interpolate( const std::vector< std::function< ValueType( const Point3D& ) > >& expr,
                     uint_t                                                             level,
                     DoFType                                                            flag = All ) const
   {
      WALBERLA_ASSERT_EQUAL( expr.size(), 1 );
      this->interpolate( expr[0], level, flag );
   };
   //@}

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
   void multElementwise( const std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >& functions,
                         uint_t                                                                             level,
                         DoFType                                                                            flag = All ) const;

   /// Replace values of the function by their inverses in an elementwise fashion
   void invertElementwise( uint_t level, DoFType flag = All, bool workOnHalos = false ) const;

   ValueType dotLocal( const EdgeDoFFunction< ValueType >& secondOp, const uint_t level, const DoFType flag = All ) const;

   /// @name Unimplemented methods (only dummys for inheritance)
   /// @{
   ValueType dotGlobal( const EdgeDoFFunction< ValueType >&, uint_t, DoFType ) const {
       WALBERLA_ABORT( "EdgeDoFFunction::dotGlobal not implemented!" )} /// @}

   ValueType sumLocal( const uint_t& level, const DoFType& flag = All, const bool& absolute = false ) const;
   ValueType sumGlobal( const uint_t& level, const DoFType& flag = All, const bool& absolute = false ) const;

   /// Set all function DoFs to zero including the ones in the halos
   void setToZero( const uint_t level ) const;

   void enumerate( uint_t level ) const;

   void enumerate( uint_t level, ValueType& offset ) const;

   const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& getVertexDataID() const { return vertexDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&   getEdgeDataID() const { return edgeDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Face >&   getFaceDataID() const { return faceDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&   getCellDataID() const { return cellDataID_; }

   ValueType getMaxValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;
   ValueType getMinValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;
   ValueType getMaxMagnitude( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;

   inline BoundaryCondition getBoundaryCondition() const { return boundaryCondition_; }
   inline void setBoundaryCondition( BoundaryCondition bc ) { boundaryCondition_ = bc; }

   template< typename OtherFunctionValueType >
   inline void copyBoundaryConditionFromFunction( const EdgeDoFFunction< OtherFunctionValueType > & other )
   {
      setBoundaryCondition( other.getBoundaryCondition() );
   }

   template < typename SenderType, typename ReceiverType >
   inline void startCommunication( const uint_t& level ) const
   {
      if ( isDummy() )
      {
         return;
      }
      communicators_.at( level )->template startCommunication< SenderType, ReceiverType >();
   }

   template < typename SenderType, typename ReceiverType >
   inline void endCommunication( const uint_t& level ) const
   {
      if ( isDummy() )
      {
         return;
      }
      communicators_.at( level )->template endCommunication< SenderType, ReceiverType >();
   }

   template < typename SenderType, typename ReceiverType >
   inline void communicate( const uint_t& level ) const
   {
      startCommunication< SenderType, ReceiverType >( level );
      endCommunication< SenderType, ReceiverType >( level );
   }

   /// Starts additive communication and will return before the communication is finished such that it can be used for
   /// asynchronous tasks. endAdditiveCommunication has to be called manually!
   /// See communicateAdditively( const uint_t& ) const for more details
   template < typename SenderType, typename ReceiverType >
   inline void startAdditiveCommunication( const uint_t& level, const bool & zeroOutDestination = true ) const
   {
      if ( isDummy() )
      {
         return;
      }
      if ( zeroOutDestination )
      {
         interpolateByPrimitiveType< ReceiverType >( real_c( 0 ), level, DoFType::All );
      }
      additiveCommunicators_.at( level )->template startCommunication< SenderType, ReceiverType >();
   }

   /// Starts additive communication that excludes primitives with a certain boundary flag from \b receiving.
   /// Will return before the communication is finished such that it can be used for
   /// asynchronous tasks. endAdditiveCommunication has to be called manually!
   /// See communicateAdditively( const uint_t&, const DoFType, const PrimitiveStorage& ) const for more details
   template < typename SenderType, typename ReceiverType >
   inline void startAdditiveCommunication( const uint_t&           level,
                                           const DoFType           boundaryTypeToSkipDuringAdditiveCommunication,
                                           const PrimitiveStorage& primitiveStorage,
                                           const bool & zeroOutDestination = true ) const
   {
      if ( isDummy() )
      {
         return;
      }
      std::vector< PrimitiveID > receiverIDs;
      std::vector< PrimitiveID > receiverNeighborIDs;
      std::vector< PrimitiveID > excludeFromReceiving;
      primitiveStorage.getPrimitiveIDsGenerically< ReceiverType >( receiverIDs );
      primitiveStorage.getNeighboringPrimitiveIDsGenerically< ReceiverType >( receiverNeighborIDs );
      for ( const PrimitiveID& id : receiverIDs )
      {
         if ( testFlag( boundaryCondition_.getBoundaryType( primitiveStorage.getPrimitive( id )->getMeshBoundaryFlag() ),
                        boundaryTypeToSkipDuringAdditiveCommunication ) )
         {
            excludeFromReceiving.push_back( id );
         }
      }
      for ( const PrimitiveID& id : receiverNeighborIDs )
      {
         if ( testFlag( boundaryCondition_.getBoundaryType( primitiveStorage.getPrimitive( id )->getMeshBoundaryFlag() ),
                        boundaryTypeToSkipDuringAdditiveCommunication ) )
         {
            excludeFromReceiving.push_back( id );
         }
      }
      if ( zeroOutDestination )
      {
         interpolateByPrimitiveType< ReceiverType >(
             real_c( 0 ), level, DoFType::All ^ boundaryTypeToSkipDuringAdditiveCommunication );
      }
      additiveCommunicators_.at( level )->template startCommunication< SenderType, ReceiverType >( excludeFromReceiving );
   }

   /// Waits for the additive communication to finish. Requires that startAdditiveCommunication() was called before.
   template < typename SenderType, typename ReceiverType >
   inline void endAdditiveCommunication( const uint_t& level ) const
   {
      if ( isDummy() )
      {
         return;
      }
      additiveCommunicators_.at( level )->template endCommunication< SenderType, ReceiverType >();
   }

   /// Additive communication sends the ghost layers of a primitive with dimension N and reduces (adds) them during the
   /// communication on the receiving primitive with dimension N-1. This is for example used for the prolongation and restriction
   /// where e.g. in 2D all the work is done on the faces and the result is additively communicated onto the edges and vertices
   /// \tparam SenderType type of the sending primitive (e.g. Face)
   /// \tparam ReceiverType type of the receiving primitive (e.g. Face)
   /// \param level the refinement level which is communicated
   /// \param zeroOutDestination if true, sets all values on the destination function to zero
   ///                           otherwise, the dst array is not modified
   template < typename SenderType, typename ReceiverType >
   inline void communicateAdditively( const uint_t& level, const bool & zeroOutDestination = true ) const
   {
      startAdditiveCommunication< SenderType, ReceiverType >( level, zeroOutDestination );
      endAdditiveCommunication< SenderType, ReceiverType >( level );
   }

   /// Similar to communicateAdditively() but excludes all primitives with the boundary type
   /// /p boundaryTypeToSkipDuringAdditiveCommunication from receiving any data. These primitives are still sending their data
   /// however!
   /// \tparam SenderType type of the sending primitive (e.g. Face)
   /// \tparam ReceiverType type of the receiving primitive (e.g. Face)
   /// \param level the refinement level which is communicated
   /// \param boundaryTypeToSkipDuringAdditiveCommunication primitives will this boundary type will no \b receive any data
   /// \param primitiveStorage
   /// \param zeroOutDestination if true, sets all values on the destination function to zero
   ///                           otherwise, the dst array is not modified
   template < typename SenderType, typename ReceiverType >
   inline void communicateAdditively( const uint_t&           level,
                                      const DoFType           boundaryTypeToSkipDuringAdditiveCommunication,
                                      const PrimitiveStorage& primitiveStorage,
                                      const bool&             zeroOutDestination = true ) const
   {
      startAdditiveCommunication< SenderType, ReceiverType >(
          level, boundaryTypeToSkipDuringAdditiveCommunication, primitiveStorage, zeroOutDestination );
      endAdditiveCommunication< SenderType, ReceiverType >( level );
   }

   /// Sets the communication mode that is used between primitives that belong to the same process. Normally this is only needed
   /// for debugging. See communication::BufferedCommunicator::LocalCommunicationMode for the available options
   /// \param localCommunicationMode
   void setLocalCommunicationMode( const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode );

   /// conversion to/from linear algebra representation
   /// @{
   void toVector( const EdgeDoFFunction< idx_t >&       numerator,
                  const std::shared_ptr< VectorProxy >& vec,
                  uint_t                                level,
                  DoFType                               flag ) const;

   void fromVector( const EdgeDoFFunction< idx_t >&       numerator,
                    const std::shared_ptr< VectorProxy >& vec,
                    uint_t                                level,
                    DoFType                               flag ) const;
   /// @}

   using Function< EdgeDoFFunction< ValueType > >::isDummy;

 private:
   template < typename PrimitiveType >
   void interpolateByPrimitiveType( const ValueType& constant, uint_t level, DoFType flag = All ) const;

   using Function< EdgeDoFFunction< ValueType > >::communicators_;
   using Function< EdgeDoFFunction< ValueType > >::additiveCommunicators_;

   PrimitiveDataID< FunctionMemory< ValueType >, Vertex > vertexDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Edge >   edgeDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Face >   faceDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Cell >   cellDataID_;

   BoundaryCondition boundaryCondition_;

   /// friend P2Function for usage of enumerate
   friend class P2Function< ValueType >;
};

// extern template class EdgeDoFFunction< double >;
extern template class EdgeDoFFunction< int >;
} // namespace hyteg
