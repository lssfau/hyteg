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

#include <memory>
#include <string>
#include <vector>

#include "core/DataTypes.h"

#include "hyteg/Function.hpp"
#include "hyteg/FunctionProperties.hpp"
#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/types/flags.hpp"
//TODO this should be improved but we need the enum which cant be forward declared
#include "hyteg/communication/BufferedCommunication.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

template < typename ValueType >
class DGFunction;
template < typename ValueType >
class FunctionMemory;
template < typename DataType, typename PrimitiveType >
class PrimitiveDataID;

class PrimitiveStorage;
class Vertex;
class Edge;
class Face;
class Cell;

namespace vertexdof {

template < typename ValueType >
class VertexDoFFunction : public Function< VertexDoFFunction< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = vertexdof::VertexDoFFunction< VType >;

   VertexDoFFunction( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage );

   VertexDoFFunction( const std::string&                         name,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      uint_t                                     minLevel,
                      uint_t                                     maxLevel );

   VertexDoFFunction( const std::string&                         name,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      uint_t                                     minLevel,
                      uint_t                                     maxLevel,
                      BoundaryCondition                          boundaryCondition );

   const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& getVertexDataID() const { return vertexDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&   getEdgeDataID() const { return edgeDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Face >&   getFaceDataID() const { return faceDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&   getCellDataID() const { return cellDataID_; }

   void swap( const VertexDoFFunction< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const;

   real_t evaluate( const Point3D& coordinates, uint_t level ) const;

   void evaluateGradient( const Point3D& coordinates, uint_t level, Point3D& gradient ) const;

   void assign( const std::vector< ValueType >&                                                      scalars,
                const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
                uint_t                                                                               level,
                DoFType                                                                              flag = All ) const;

   void assign( const P2Function< ValueType >& src, const uint_t& P1Level, const DoFType& flag = All ) const;

   void add( const ValueType& scalar, const uint_t& level, DoFType flag = All ) const;

   void add( const std::vector< ValueType >&                                                      scalars,
             const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
             uint_t                                                                               level,
             DoFType                                                                              flag = All ) const;

   void multElementwise( const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
                         uint_t                                                                               level,
                         DoFType                                                                              flag = All ) const;

   ValueType dotLocal( const VertexDoFFunction< ValueType >& rhs, uint_t level, DoFType flag = All ) const;
   ValueType dotGlobal( const VertexDoFFunction< ValueType >& rhs, uint_t level, DoFType flag = All ) const;

   ValueType sumLocal( const uint_t& level, const DoFType& flag = All, const bool & absolute = false ) const;
   ValueType sumGlobal( const uint_t& level, const DoFType& flag = All, const bool & absolute = false ) const;

   void integrateDG( DGFunction< ValueType >& rhs, VertexDoFFunction< ValueType >& rhsP1, uint_t level, DoFType flag );

   /// Interpolates a given expression to a VertexDoFFunction
   void interpolate( const ValueType& constant, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, BoundaryUID boundaryUID ) const;

   void interpolateExtended( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                             const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& srcFunctions,
                             uint_t                                                                               level,
                             DoFType flag = All ) const;

   void interpolateExtended( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                             const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& srcFunctions,
                             uint_t                                                                               level,
                             BoundaryUID boundaryUID ) const;

   /// assigns unique values to all data points
   /// this function is mainly used for petsc to get global identifier for all DoFs
   /// \tparam ValueType
   /// \param level
   void enumerate( uint_t level ) const;

   // TODO: write more general version(s)
   ValueType getMaxValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;
   ValueType getMinValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;
   ValueType getMaxMagnitude( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;

   BoundaryCondition getBoundaryCondition() const;

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
   inline void startAdditiveCommunication( const uint_t& level ) const
   {
      if ( isDummy() )
      {
         return;
      }
      interpolateByPrimitiveType< ReceiverType >(
          real_c( 0 ), level, DoFType::All );
      additiveCommunicators_.at( level )->template startCommunication< SenderType, ReceiverType >();
   }

   /// Starts additive communication that excludes primitives with a certain boundary flag from \b receiving.
   /// Will return before the communication is finished such that it can be used for
   /// asynchronous tasks. endAdditiveCommunication has to be called manually!
   /// See communicateAdditively( const uint_t&, const DoFType, const PrimitiveStorage& ) const for more details
   template < typename SenderType, typename ReceiverType >
   inline void startAdditiveCommunication( const uint_t&           level,
                                           const DoFType           boundaryTypeToSkipDuringAdditiveCommunication,
                                           const PrimitiveStorage& primitiveStorage ) const
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
      interpolateByPrimitiveType< ReceiverType >(
          real_c( 0 ), level, DoFType::All ^ boundaryTypeToSkipDuringAdditiveCommunication_ );
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
   template < typename SenderType, typename ReceiverType >
   inline void communicateAdditively( const uint_t& level ) const
   {
      startAdditiveCommunication< SenderType, ReceiverType >( level );
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
   template < typename SenderType, typename ReceiverType >
   inline void communicateAdditively( const uint_t&           level,
                                      const DoFType           boundaryTypeToSkipDuringAdditiveCommunication,
                                      const PrimitiveStorage& primitiveStorage ) const
   {
      startAdditiveCommunication< SenderType, ReceiverType >( level, boundaryTypeToSkipDuringAdditiveCommunication, primitiveStorage );
      endAdditiveCommunication< SenderType, ReceiverType >( level );
   }

   /// Sets the communication mode that is used between primitives that belong to the same process. Normally this is only needed
   /// for debugging. See communication::BufferedCommunicator::LocalCommunicationMode for the available options
   /// \param localCommunicationMode
   void setLocalCommunicationMode( const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode );

   using Function< VertexDoFFunction< ValueType > >::isDummy;

 private:
   template < typename PrimitiveType >
   void interpolateByPrimitiveType( const ValueType& constant, uint_t level, DoFType flag = All ) const;

   void enumerate( uint_t level, ValueType& offset ) const;

   using Function< VertexDoFFunction< ValueType > >::communicators_;
   using Function< VertexDoFFunction< ValueType > >::additiveCommunicators_;

   PrimitiveDataID< FunctionMemory< ValueType >, Vertex > vertexDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Edge >   edgeDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Face >   faceDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Cell >   cellDataID_;

   BoundaryCondition boundaryCondition_;

   DoFType boundaryTypeToSkipDuringAdditiveCommunication_;

   /// friend Stokes and P2Function for usage of enumerate
   friend class P2Function< ValueType >;
   friend class P1StokesFunction< ValueType >;
   friend class P2P1TaylorHoodFunction< ValueType >;
};

inline void projectMean( const VertexDoFFunction< real_t >& pressure, const uint_t& level )
{
   if ( pressure.isDummy() )
   {
      return;
   }
   const uint_t numGlobalVertices = numberOfGlobalDoFs< VertexDoFFunctionTag >( *pressure.getStorage(), level );
   const real_t sum               = pressure.sumGlobal( level, All );
   pressure.add( -sum / real_c( numGlobalVertices ), level, All );
}

} // namespace vertexdof
} // namespace hyteg
