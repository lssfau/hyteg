/*
 * Copyright (c) 2017-2021 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"

#include "DGDataHandling.hpp"
#include "DGEdge.hpp"
#include "DGFace.hpp"
#include "DGMemory.hpp"
#include "DGPackInfo.hpp"
#include "DGVertex.hpp"

namespace hyteg {

template < typename ValueType >
class DGFunction final : public Function< DGFunction< ValueType > >
{
 public:
   DGFunction( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : DGFunction( name, storage, minLevel, maxLevel, BoundaryCondition::create0123BC() )
   {}

   DGFunction( const std::string&                         name,
               const std::shared_ptr< PrimitiveStorage >& storage,
               uint_t                                     minLevel,
               uint_t                                     maxLevel,
               BoundaryCondition                          boundaryCondition )
   : Function< DGFunction< ValueType > >( name, storage, minLevel, maxLevel )
   , boundaryCondition_( boundaryCondition )
   {
      auto vertexDGFunctionMemoryDataHandling =
          std::make_shared< VertexDGFunctionMemoryDataHandling< ValueType > >( minLevel, maxLevel );
      auto edgeDGFunctionMemoryDataHandling =
          std::make_shared< EdgeDGFunctionMemoryDataHandling< ValueType > >( minLevel, maxLevel );
      auto faceDGFunctionMemoryDataHandling =
          std::make_shared< FaceDGFunctionMemoryDataHandling< ValueType > >( minLevel, maxLevel );

      storage->addFaceData( faceDataID_, faceDGFunctionMemoryDataHandling, name );
      storage->addEdgeData( edgeDataID_, edgeDGFunctionMemoryDataHandling, name );
      storage->addVertexData( vertexDataID_, vertexDGFunctionMemoryDataHandling, name );
      for ( uint_t level = minLevel; level <= maxLevel; ++level )
      {
         //communicators_[level]->setLocalCommunicationMode(communication::BufferedCommunicator::BUFFERED_MPI);
         communicators_[level]->addPackInfo(
             std::make_shared< DGPackInfo< ValueType > >( level, vertexDataID_, edgeDataID_, faceDataID_, this->getStorage() ) );
      }
   }

   inline void interpolate( std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType flag = All );

   inline void interpolate( const ValueType& constant, uint_t level, DoFType flag = All );

   inline void interpolateExtended( std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                                    const std::vector< DGFunction< ValueType >* >                                  srcFunctions,
                                    uint_t                                                                         level,
                                    DoFType                                                                        flag = All );

   inline void assign( const std::vector< ValueType >&                                               scalars,
                       const std::vector< std::reference_wrapper< const DGFunction< ValueType > > >& functions,
                       uint_t                                                                        level,
                       DoFType                                                                       flag = All );

   inline void add( const std::vector< ValueType >                scalars,
                    const std::vector< DGFunction< ValueType >* > functions,
                    uint_t                                        level,
                    DoFType                                       flag = All );
   inline void enumerate( uint_t level, ValueType offset );

   inline void enumerate( uint_t level );

   inline real_t getMaxValue( const uint_t level, DoFType flag = All );
   inline real_t getMinValue( const uint_t level, DoFType flag = All );
   inline real_t getMaxMagnitude( const uint_t level, DoFType flag = All );

   const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& getVertexDataID() const { return vertexDataID_; }

   const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& getEdgeDataID() const { return edgeDataID_; }

   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& getFaceDataID() const { return faceDataID_; }

   void projectP1( P1Function< real_t >& src, uint_t level, DoFType flag, UpdateType updateType = Replace );

   inline BoundaryCondition getBoundaryCondition() const { return boundaryCondition_; }

   template < typename SenderType, typename ReceiverType >
   inline void startCommunication( const uint_t& level ) const
   {
      communicators_.at( level )->template startCommunication< SenderType, ReceiverType >();
   }

   template < typename SenderType, typename ReceiverType >
   inline void endCommunication( const uint_t& level ) const
   {
      communicators_.at( level )->template endCommunication< SenderType, ReceiverType >();
   }

   template < typename SenderType, typename ReceiverType >
   inline void communicate( const uint_t& level ) const
   {
      startCommunication< SenderType, ReceiverType >( level );
      endCommunication< SenderType, ReceiverType >( level );
   }

   inline void
       setLocalCommunicationMode( const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode )
   {
      for ( auto& communicator : communicators_ )
      {
         communicator.second->setLocalCommunicationMode( localCommunicationMode );
      }
   }

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
   inline void multElementwise( const std::vector< std::reference_wrapper< const DGFunction< ValueType > > >& functions,
                                uint_t                                                                        level,
                                DoFType                                                                       flag = All ) const;

   void add( const ValueType scalar, uint_t level, DoFType flag = All ) const;

   void copyFrom( const DGFunction< ValueType >& other, const uint_t& level ) const
   {
      WALBERLA_ABORT( "DGFunction::copyFrom not implemented!" )
   }

   void copyFrom( const DGFunction< ValueType >&            other,
                  const uint_t&                                  level,
                  const std::map< PrimitiveID::IDType, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID::IDType, uint_t >& otherPrimitiveIDsToRank ) const
   {
      WALBERLA_ABORT( "DGFunction::copyFrom not implemented!" )
   }

   void interpolate( ValueType constant, uint_t level, DoFType flag = All ) const
   {
      WALBERLA_ABORT( "DGFunction::interpolate not implemented!" )
   }

   void interpolate( const std::function< ValueType( const hyteg::Point3D& ) >& expr, uint_t level, DoFType flag = All ) const
   {
      WALBERLA_ABORT( "DGFunction::interpolate not implemented!" )
   }

   void interpolate( const std::vector< std::function< ValueType( const hyteg::Point3D& ) > >& expr,
                     uint_t                                                                    level,
                     DoFType                                                                   flag = All ) const
   {
      WALBERLA_ABORT( "DGFunction::interpolate not implemented!" )
   }

   void setBoundaryCondition( BoundaryCondition bc ){WALBERLA_ABORT( "DGFunction::setBoundaryCondition not implemented!" )}

   ValueType dotLocal( const DGFunction< ValueType >& secondOp, uint_t level, DoFType flag ) const
   {
      WALBERLA_ABORT( "DGFunction::dotLocal not implemented!" )
   }

   ValueType dotGlobal( const DGFunction< ValueType >& secondOp, uint_t level, DoFType flag ) const
   {
      WALBERLA_ABORT( "DGFunction::dotGlobal not implemented!" )
   }

   void assign( const std::vector< ValueType >&                                               scalars,
                const std::vector< std::reference_wrapper< const DGFunction< ValueType > > >& functions,
                uint_t                                                                        level,
                DoFType                                                                       flag = All ) const
   {
      WALBERLA_ABORT( "DGFunction::assign with const qualifier not implemented!" )
   }

   void swap( const DGFunction< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const
   {
      WALBERLA_ABORT( "DGFunction::swap not implemented!" )
   }

   void add( const std::vector< ValueType >&                                                    scalars,
             const std::vector< std::reference_wrapper< const DGFunction< ValueType > > >& functions,
             uint_t                                                                             level,
             DoFType                                                                            flag = All ) const
   {
      WALBERLA_ABORT( "DGFunction::add for vector of functions not implemented!" )
   }

   /// @}

 private:
   using Function< DGFunction< ValueType > >::communicators_;

   PrimitiveDataID< FunctionMemory< ValueType >, Vertex > vertexDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Edge >   edgeDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Face >   faceDataID_;

   BoundaryCondition boundaryCondition_;
};

template < typename ValueType >
void DGFunction< ValueType >::add( const std::vector< ValueType >                scalars,
                                   const std::vector< DGFunction< ValueType >* > functions,
                                   uint_t                                        level,
                                   DoFType                                       flag )
{
   // Collect all source IDs in a vector
   //  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Vertex>> srcVertexIDs;
   //  std::vector<PrimitiveDataID<FunctionMemory< ValueType >, Edge>>   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > srcFaceIDs;

   for ( auto& function : functions )
   {
      //    srcVertexIDs.push_back(function->vertexDataID_);
      //    srcEdgeIDs.push_back(function->edgeDataID_);
      srcFaceIDs.push_back( function->faceDataID_ );
   }

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         DGFace::add< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_ );
      }
   }
}

template < typename ValueType >
inline void DGFunction< ValueType >::interpolate( std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType flag )
{
   std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) > exprExtended =
       [&expr]( const hyteg::Point3D& x, const std::vector< ValueType >& ) { return expr( x ); };
   interpolateExtended( exprExtended, {}, level, flag );
}

template < typename ValueType >
inline void DGFunction< ValueType >::interpolate( const ValueType& constant, uint_t level, DoFType flag )
{
   std::function< ValueType( const Point3D& ) > auxFunc = [constant]( const hyteg::Point3D& x ) { return constant; };
   this->interpolate( {auxFunc}, level, flag );
}

template < typename ValueType >
void DGFunction< ValueType >::interpolateExtended(
    std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
    const std::vector< DGFunction< ValueType >* >                                  srcFunctions,
    uint_t                                                                         level,
    DoFType                                                                        flag )
{
   this->startTiming( "Interpolate" );
   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;

   for ( auto& function : srcFunctions )
   {
      srcVertexIDs.push_back( function->vertexDataID_ );
      srcEdgeIDs.push_back( function->edgeDataID_ );
      srcFaceIDs.push_back( function->faceDataID_ );
   }

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         DGVertex::interpolate< ValueType >( level, vertex, vertexDataID_, srcVertexIDs, expr, this->getStorage() );
      }
   }

   communicators_[level]->template startCommunication< Vertex, Edge >();

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         DGEdge::interpolate< ValueType >( level, edge, edgeDataID_, srcEdgeIDs, expr, this->getStorage() );
      }
   }

   communicators_[level]->template endCommunication< Vertex, Edge >();
   communicators_[level]->template startCommunication< Edge, Face >();

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         DGFace::interpolate< ValueType >( level, face, faceDataID_, srcFaceIDs, expr );
      }
   }

   communicators_[level]->template endCommunication< Edge, Face >();
   this->stopTiming( "Interpolate" );
}

template < typename ValueType >
void DGFunction< ValueType >::assign( const std::vector< ValueType >&                                               scalars,
                                      const std::vector< std::reference_wrapper< const DGFunction< ValueType > > >& functions,
                                      uint_t                                                                        level,
                                      DoFType                                                                       flag )
{
   this->startTiming( "Assign" );

   WALBERLA_ASSERT_EQUAL( scalars.size(), functions.size() )

   // Collect all source IDs in a vector
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcVertexIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > >   srcEdgeIDs;
   std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > >   srcFaceIDs;

   for ( const DGFunction< ValueType >& function : functions )
   {
      srcVertexIDs.push_back( function.vertexDataID_ );
      srcEdgeIDs.push_back( function.edgeDataID_ );
      srcFaceIDs.push_back( function.faceDataID_ );
   }

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         DGVertex::assign< ValueType >( level, vertex, scalars, srcVertexIDs, vertexDataID_ );
      }
   }

   communicators_[level]->template startCommunication< Vertex, Edge >();

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         DGEdge::assign< ValueType >( level, edge, scalars, srcEdgeIDs, edgeDataID_ );
      }
   }

   communicators_[level]->template endCommunication< Vertex, Edge >();
   communicators_[level]->template startCommunication< Edge, Face >();

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         DGFace::assign< ValueType >( level, face, scalars, srcFaceIDs, faceDataID_ );
      }
   }

   communicators_[level]->template endCommunication< Edge, Face >();
   this->stopTiming( "Assign" );
}

template < typename ValueType >
void DGFunction< ValueType >::enumerate( uint_t level )
{
   enumerate( level, static_cast< ValueType >( 0 ) );
}

template < typename ValueType >
void DGFunction< ValueType >::enumerate( uint_t level, ValueType offset )
{
   this->startTiming( "Enumerate" );

   uint_t counter = hyteg::numberOfLocalDoFs< VertexDoFFunctionTag >( *( this->getStorage() ), level );

   std::vector< uint_t > dofs_per_rank = walberla::mpi::allGather( counter );

   // the next line does not make sense, because the implementation of enumerate()
   // in DGVertex, DGEdge and DGFace expects to receive a uint_t!
   // ValueType startOnRank = offset;
   // replaced it by this "hotfix"
   auto startOnRank = static_cast< uint_t >( offset );

   for ( uint_t i = 0; i < uint_c( walberla::MPIManager::instance()->rank() ); ++i )
   {
      startOnRank += dofs_per_rank[i];
   }
   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;
      DGVertex::enumerate( vertex, vertexDataID_, level, startOnRank );
   }

   communicators_[level]->template startCommunication< Vertex, Edge >();
   communicators_[level]->template endCommunication< Vertex, Edge >();

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;
      DGEdge::enumerate< ValueType >( level, edge, edgeDataID_, startOnRank );
   }

   communicators_[level]->template startCommunication< Edge, Face >();
   communicators_[level]->template endCommunication< Edge, Face >();

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;
      DGFace::enumerate< ValueType >( level, face, faceDataID_, startOnRank );
   }

   communicators_[level]->template startCommunication< Face, Edge >();
   communicators_[level]->template endCommunication< Face, Edge >();

   communicators_[level]->template startCommunication< Edge, Vertex >();
   communicators_[level]->template endCommunication< Edge, Vertex >();
   this->stopTiming( "Enumerate" );
}

template < typename ValueType >
void DGFunction< ValueType >::projectP1( P1Function< real_t >& src, uint_t level, DoFType flag, UpdateType updateType )
{
   this->startTiming( "projectP1" );

   src.startCommunication< Edge, Vertex >( level );
   src.startCommunication< Face, Edge >( level );
   src.endCommunication< Edge, Vertex >( level );

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         DGVertex::projectP1< real_t >(
             level, vertex, this->getStorage(), src.getVertexDataID(), this->getVertexDataID(), updateType );
      }
   }

   startCommunication< Vertex, Edge >( level );

   src.endCommunication< Face, Edge >( level );

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );

      if ( testFlag( edgeBC, flag ) )
      {
         DGEdge::projectP1< real_t >( level, edge, this->getStorage(), src.getEdgeDataID(), this->getEdgeDataID(), updateType );
      }
   }

   endCommunication< Vertex, Edge >( level );

   startCommunication< Edge, Face >( level );

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         DGFace::projectP1< real_t >( level, face, this->getStorage(), src.getFaceDataID(), this->getFaceDataID(), updateType );
      }
   }

   endCommunication< Edge, Face >( level );

   this->stopTiming( "projectP1" );
}

template < typename ValueType >
real_t DGFunction< ValueType >::getMaxValue( const uint_t level, DoFType flag )
{
   real_t localMax = -std::numeric_limits< ValueType >::max();

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         localMax = std::max( localMax, DGFace::getMaxValue< ValueType >( level, face, faceDataID_ ) );
      }
   }

   walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );
   return localMax;
}

template < typename ValueType >
real_t DGFunction< ValueType >::getMinValue( const uint_t level, DoFType flag )
{
   ValueType localMin = std::numeric_limits< ValueType >::max();

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         localMin = std::min( localMin, DGFace::getMinValue< ValueType >( level, face, faceDataID_ ) );
      }
   }

   walberla::mpi::allReduceInplace( localMin, walberla::mpi::MIN, walberla::mpi::MPIManager::instance()->comm() );
   return localMin;
}

template < typename ValueType >
real_t DGFunction< ValueType >::getMaxMagnitude( const uint_t level, DoFType flag )
{
   ValueType localMax = -std::numeric_limits< ValueType >::max();

   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face&         face   = *it.second;
      const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         localMax = std::max( localMax, DGFace::getMaxMagnitude< ValueType >( level, face, faceDataID_ ) );
      }
   }

   walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );
   return localMax;
}

template < typename ValueType >
void DGFunction< ValueType >::multElementwise(
    const std::vector< std::reference_wrapper< const DGFunction< ValueType > > >& functions,
    uint_t                                                                        level,
    DoFType                                                                       flag ) const
{
   this->startTiming( "Multiply elementwise" );

   if ( this->getStorage()->hasGlobalCells() )
   {
      WALBERLA_ABORT( "DGFunction::multElementwise() not implemented for 3D!" );
   }

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Vertex > > srcIDs;
      for ( const DGFunction& function : functions )
      {
         srcIDs.push_back( function.getVertexDataID() );
      }

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         DGVertex::multElementwise( level, vertex, srcIDs, vertexDataID_ );
      }
   }
   for ( auto& it : this->getStorage()->getEdges() )
   {
      std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Edge > > srcIDs;
      for ( const DGFunction& function : functions )
      {
         srcIDs.push_back( function.getEdgeDataID() );
      }

      Edge& edge = *it.second;
      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         DGEdge::multElementwise( level, edge, srcIDs, edgeDataID_ );
      }
   }
   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      std::vector< PrimitiveDataID< FunctionMemory< ValueType >, Face > > srcIDs;
      for ( const DGFunction& function : functions )
      {
         srcIDs.push_back( function.getFaceDataID() );
      }
      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         DGFace::multElementwise( level, face, srcIDs, faceDataID_ );
      }
   }
   this->stopTiming( "Multiply elementwise" );
}

template < typename ValueType >
void DGFunction< ValueType >::add( const ValueType scalar, uint_t level, DoFType flag ) const
{
   this->startTiming( "Add scalar" );

   if ( this->getStorage()->hasGlobalCells() )
   {
      WALBERLA_ABORT( "DGFunction::add() not implemented for 3D!" );
   }

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
      {
         DGVertex::add( level, vertex, scalar, vertexDataID_ );
      }
   }
   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
      {
         DGEdge::add( level, edge, scalar, edgeDataID_ );
      }
   }
   for ( auto& it : this->getStorage()->getFaces() )
   {
      Face& face = *it.second;

      if ( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
      {
         DGFace::add( level, face, scalar, faceDataID_ );
      }
   }
   this->stopTiming( "Add scalar" );
}

} // namespace hyteg
