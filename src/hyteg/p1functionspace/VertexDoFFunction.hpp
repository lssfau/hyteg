/*
 * Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl, Benjamin Mann.
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

#include "hyteg/ReferenceCounter.hpp"
#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/memory/FunctionMemory.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"
#include "hyteg/types/types.hpp"
/// \todo This should be improved, but we need the enum which can't be forward declared
#include "hyteg/communication/BufferedCommunication.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

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
class VertexDoFFunction final : public Function< VertexDoFFunction< ValueType > >
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



    VertexDoFFunction( const std::string&                         name,
                       const std::shared_ptr< PrimitiveStorage >& storage,
                       uint_t                                     minLevel,
                       uint_t                                     maxLevel,
                       BoundaryCondition                          boundaryCondition,
                       bool                                       addVolumeGhostLayer );

   ~VertexDoFFunction();

   /// Copy constructor
   VertexDoFFunction( const VertexDoFFunction< ValueType >& other );

   /// Copy assignment
   VertexDoFFunction& operator=( const VertexDoFFunction< ValueType >& other );

   virtual uint_t getDimension() const { return 1; }

   bool hasMemoryAllocated( const uint_t& level, const Vertex& vertex ) const;
   bool hasMemoryAllocated( const uint_t& level, const Edge& edge ) const;
   bool hasMemoryAllocated( const uint_t& level, const Face& face ) const;
   bool hasMemoryAllocated( const uint_t& level, const Cell& cell ) const;


   void allocateMemory( const uint_t& level, const Vertex& vertex );
   void allocateMemory( const uint_t& level, const Edge& edge );
   void allocateMemory( const uint_t& level, const Face& face );
   void allocateMemory( const uint_t& level, const Cell& cell );

   void deleteMemory( const uint_t& level, const Vertex& vertex );
   void deleteMemory( const uint_t& level, const Edge& edge );
   void deleteMemory( const uint_t& level, const Face& face );
   void deleteMemory( const uint_t& level, const Cell& cell );

   const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& getVertexDataID() const { return vertexDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&   getEdgeDataID() const { return edgeDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Face >&   getFaceDataID() const { return faceDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&   getCellDataID() const { return cellDataID_; }

   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& getFaceGLDataID( uint_t glID ) const
   {
      WALBERLA_CHECK( hasVolumeGhostLayer() );
      WALBERLA_CHECK( !this->getStorage()->hasGlobalCells() );
      return faceGhostLayerDataIDs_.at( glID );
   }
   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >& getCellGLDataID( uint_t glID ) const
   {
      WALBERLA_CHECK( hasVolumeGhostLayer() );
      WALBERLA_CHECK( this->getStorage()->hasGlobalCells() );
      return cellGhostLayerDataIDs_.at( glID );
   }

   void swap( const VertexDoFFunction< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const;

   /// \brief Copies all values function data from other to this.
   ///
   /// This method can be used safely if the other function is located on a different PrimitiveStorage.
   /// Both storages must have identical distribution.
   ///
   template< typename otherValueType >
   void copyFrom( const VertexDoFFunction< otherValueType >& other, const uint_t& level ) const;

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
   void copyFrom( const VertexDoFFunction< ValueType >&  other,
                  const uint_t&                          level,
                  const std::map< PrimitiveID, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID, uint_t >& otherPrimitiveIDsToRank ) const;

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
   /// \param id PrimitiveID of the cell/face where physicalCoords is located. If a valid id is provided, the loop over all local primitives is skipped.
   /// \return true if the function was evaluated successfully, false otherwise
   ///
   bool evaluate( const Point3D& physicalCoords,
                  uint_t         level,
                  ValueType&     value,
                  real_t         searchToleranceRadius = real_c(1e-05),
                  PrimitiveID    id                    = PrimitiveID() ) const;

   void evaluateGradient( const Point3D& physicalCoords, uint_t level, Point3D& gradient ) const;

   void assign( const std::vector< ValueType >&                                                      scalars,
                const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
                uint_t                                                                               level,
                DoFType                                                                              flag = All ) const;

   void add( ValueType scalar, uint_t level, DoFType flag = All ) const;

   void add( const std::vector< ValueType >&                                                      scalars,
             const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
             uint_t                                                                               level,
             DoFType                                                                              flag = All ) const;

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
   void multElementwise( const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
                         uint_t                                                                               level,
                         DoFType                                                                              flag = All ) const;

   /// Replace values of the function by their inverses in an elementwise fashion
   void invertElementwise( uint_t level, DoFType flag = All, bool workOnHalos = false ) const;

   ValueType dotLocal( const VertexDoFFunction< ValueType >& rhs, uint_t level, DoFType flag = All ) const;
   ValueType dotGlobal( const VertexDoFFunction< ValueType >& rhs, uint_t level, DoFType flag = All ) const;

   ValueType sumLocal( const uint_t& level, const DoFType& flag = All, const bool& absolute = false ) const;
   ValueType sumGlobal( const uint_t& level, const DoFType& flag = All, const bool& absolute = false ) const;

   /// @name Member functions for interpolation using BoundaryUID flags
   //@{
   void interpolate( ValueType constant, uint_t level, BoundaryUID boundaryUID ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, BoundaryUID boundaryUID ) const;

   void interpolate( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                     const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& srcFunctions,
                     uint_t                                                                               level,
                     BoundaryUID                                                                          boundaryUID ) const;
   //@}

   /// @name Member functions for interpolation using DoFType flags
   //@{
   /// Interpolates a given expression to a VertexDoFFunction
   void interpolate( ValueType constant, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::vector< std::function< ValueType( const Point3D& ) > >& expr,
                     uint_t                                                             level,
                     DoFType                                                            flag = All ) const
   {
      WALBERLA_ASSERT_EQUAL( expr.size(), 1 );
      this->interpolate( expr[0], level, flag );
   };

   void interpolate( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                     const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& srcFunctions,
                     uint_t                                                                               level,
                     DoFType                                                                              flag = All ) const;
   //@}

   /// interpolate data from a coarser mesh.
   ///
   /// Note that this only works for adaptively refined meshes and only if the
   /// child-primitives are located on the same processes as their parents!
   /// Involves global communication!
   void interpolate( const VertexDoFFunction< ValueType >& functionOnParentGrid, uint_t level, uint_t srcLevel );

   /// Set all function DoFs to zero including the ones in the halos
   void setToZero( const uint_t level ) const;

   /// assigns unique values to all data points
   /// this function is mainly used for petsc to get global identifier for all DoFs
   /// \tparam ValueType
   /// \param level
   void enumerate( uint_t level ) const;

   /// like enumerate() but starting with value given by offset parameter
   void enumerate( uint_t level, ValueType& offset ) const;

   // TODO: write more general version(s)
   ValueType getMaxValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;
   ValueType getMinValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;
   ValueType getMaxMagnitude( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;

   BoundaryCondition getBoundaryCondition() const;
   void              setBoundaryCondition( BoundaryCondition bc );

   template < typename OtherFunctionValueType >
   void copyBoundaryConditionFromFunction( const VertexDoFFunction< OtherFunctionValueType >& other )
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
   inline void startAdditiveCommunication( const uint_t& level, const bool& zeroOutDestination = true ) const
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
                                           const bool&             zeroOutDestination = true ) const
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
   inline void communicateAdditively( const uint_t& level, const bool& zeroOutDestination = true ) const
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
   void toVector( const VertexDoFFunction< idx_t >&     numerator,
                  const std::shared_ptr< VectorProxy >& vec,
                  uint_t                                level,
                  DoFType                               flag ) const;

   void fromVector( const VertexDoFFunction< idx_t >&     numerator,
                    const std::shared_ptr< VectorProxy >& vec,
                    uint_t                                level,
                    DoFType                               flag ) const;
   /// @}

   using Function< VertexDoFFunction< ValueType > >::isDummy;

   /// \brief Returns the number of DoFs that are allocated on this process.
   uint_t getNumberOfLocalDoFs( uint_t level ) const { return numberOfLocalDoFs< P1FunctionTag >( *this->storage_, level ); }

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
      return numberOfGlobalDoFs< P1FunctionTag >( *this->storage_, level, communicator, onRootOnly );
   }

   bool hasVolumeGhostLayer() const
   {
      if ( this->getStorage()->hasGlobalCells() )
      {
         return !cellGhostLayerDataIDs_.empty();
      }
      else
      {
         return !faceGhostLayerDataIDs_.empty();
      }
   }


 private:
   template < typename PrimitiveType >
   void interpolateByPrimitiveType( const ValueType& constant, uint_t level, DoFType flag = All ) const;

   inline void deleteFunctionMemory()
   {
      this->storage_->deleteVertexData( vertexDataID_ );
      this->storage_->deleteEdgeData( edgeDataID_ );
      this->storage_->deleteFaceData( faceDataID_ );
      this->storage_->deleteCellData( cellDataID_ );
   }

   using Function< VertexDoFFunction< ValueType > >::communicators_;
   using Function< VertexDoFFunction< ValueType > >::additiveCommunicators_;

   PrimitiveDataID< FunctionMemory< ValueType >, Vertex > vertexDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Edge >   edgeDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Face >   faceDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Cell >   cellDataID_;

   BoundaryCondition boundaryCondition_;

   /// All functions that actually allocate data and are not composites are handles to the allocated memory.
   /// This means the copy-ctor and copy-assignment only create a handle that is associated with the same memory.
   /// Deep copies must be created explicitly.
   /// To make sure that functions that are not used anymore are deleted, we need to add this reference counter to the handle.
   /// Once it drops to zero, we can deallocate the memory from the storage.
   std::shared_ptr< internal::ReferenceCounter > referenceCounter_;

   /// One data ID per ghost-layer. Should be up to 3 in 2D and 4 in 3D.
   /// Maps from the local macro-edge ID (for 2D) or local macro-face ID (for 3D) to the respective ghost-layer memory.
   std::map< uint_t, PrimitiveDataID< FunctionMemory< ValueType >, Face > > faceGhostLayerDataIDs_;
   std::map< uint_t, PrimitiveDataID< FunctionMemory< ValueType >, Cell > > cellGhostLayerDataIDs_;
};

inline void projectMean( const VertexDoFFunction< real_t >& pressure, const uint_t& level )
{
   if ( pressure.isDummy() )
   {
      return;
   }
   const uint_t numGlobalVertices = numberOfGlobalDoFs< VertexDoFFunctionTag >(
       *pressure.getStorage(), level, pressure.getStorage()->getSplitCommunicatorByPrimitiveDistribution() );
   const real_t sum = pressure.sumGlobal( level, All );
   pressure.add( -sum / real_c( numGlobalVertices ), level, All );
}

// extern template class VertexDoFFunction< double >;
extern template class VertexDoFFunction< int >;
extern template class VertexDoFFunction< idx_t >;

template < typename ValueType >
template < typename otherValueType >
void VertexDoFFunction< ValueType >::copyFrom( const VertexDoFFunction< otherValueType >& other, const uint_t& level ) const
{
   if ( isDummy() )
   {
      return;
   }
   this->startTiming( "Copy" );

   for ( auto& it : this->getStorage()->getVertices() )
   {
      auto primitiveID = it.first;
      WALBERLA_ASSERT( other.getStorage()->vertexExistsLocally( primitiveID ) )
      this->getStorage()
          ->getVertex( primitiveID )
          ->getData( vertexDataID_ )
          ->copyFrom( *other.getStorage()->getVertex( primitiveID )->getData( other.getVertexDataID() ), level );
   }

   for ( auto& it : this->getStorage()->getEdges() )
   {
      auto primitiveID = it.first;
      WALBERLA_ASSERT( other.getStorage()->edgeExistsLocally( primitiveID ) )
      this->getStorage()
          ->getEdge( primitiveID )
          ->getData( edgeDataID_ )
          ->copyFrom( *other.getStorage()->getEdge( primitiveID )->getData( other.getEdgeDataID() ), level );
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      auto primitiveID = it.first;
      WALBERLA_ASSERT( other.getStorage()->faceExistsLocally( primitiveID ) )
      this->getStorage()
          ->getFace( primitiveID )
          ->getData( faceDataID_ )
          ->copyFrom( *other.getStorage()->getFace( primitiveID )->getData( other.getFaceDataID() ), level );
   }

   for ( auto& it : this->getStorage()->getCells() )
   {
      auto primitiveID = it.first;
      WALBERLA_ASSERT( other.getStorage()->cellExistsLocally( primitiveID ) )
      this->getStorage()
          ->getCell( primitiveID )
          ->getData( cellDataID_ )
          ->copyFrom( *other.getStorage()->getCell( primitiveID )->getData( other.getCellDataID() ), level );
   }

   this->stopTiming( "Copy" );
}

} // namespace vertexdof
} // namespace hyteg