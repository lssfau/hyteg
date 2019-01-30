#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFFunction.hpp"
#include "tinyhhg_core/p2functionspace/P2Multigrid.hpp"
#include "tinyhhg_core/p2functionspace/P2TransferOperators.hpp"

namespace hhg {

using walberla::real_c;

template < typename ValueType >
class P2Function : public Function< P2Function< ValueType > >
{
 public:

  typedef ValueType valueType;

  template< typename VType >
  using FunctionType = P2Function< VType >;

   P2Function( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage )
   : Function< P2Function< ValueType > >( name, storage )
   , vertexDoFFunction_( vertexdof::VertexDoFFunction< ValueType >( name + "_VertexDoF_dummy", storage ) )
   , edgeDoFFunction_( EdgeDoFFunction< ValueType >( name + "__EdgeDoF_dummy", storage ) )
   {}

   P2Function( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : P2Function( name, storage, minLevel, maxLevel, BoundaryCondition::create012BC() )
   {}

   P2Function( const std::string&                         name,
               const std::shared_ptr< PrimitiveStorage >& storage,
               uint_t                                     minLevel,
               uint_t                                     maxLevel,
               BoundaryCondition                          boundaryCondition )
   : Function< P2Function< ValueType > >( name, storage, minLevel, maxLevel )
   , vertexDoFFunction_(
         vertexdof::VertexDoFFunction< ValueType >( name + "_VertexDoF", storage, minLevel, maxLevel, boundaryCondition ) )
   , edgeDoFFunction_( EdgeDoFFunction< ValueType >( name + "_EdgeDoF", storage, minLevel, maxLevel, boundaryCondition ) )
   {
      for( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         /// one has to use the communicators of the vertexDoF and edgeDoF function to communicate
         /// TODO: find better solution
         communicators_[level] = NULL;
      }
   }

   vertexdof::VertexDoFFunction< ValueType > getVertexDoFFunctionCopy() const { return vertexDoFFunction_; }
   EdgeDoFFunction< ValueType >              getEdgeDoFFunctionCopy() const { return edgeDoFFunction_; }

   const vertexdof::VertexDoFFunction< ValueType > & getVertexDoFFunction() const { return vertexDoFFunction_; }
   const EdgeDoFFunction< ValueType > &              getEdgeDoFFunction() const { return edgeDoFFunction_; }

    template < typename SenderType, typename ReceiverType >
    inline void communicate( const uint_t& level ) const
    {
      vertexDoFFunction_.template communicate< SenderType, ReceiverType >( level );
      edgeDoFFunction_  .template communicate< SenderType, ReceiverType >( level );
    }

    inline void interpolate( const ValueType& constant, uint_t level, DoFType flag = All ) const
   {
      vertexDoFFunction_.interpolate( constant, level, flag );
      edgeDoFFunction_.interpolate( constant, level, flag );
   }

   inline void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType flag = All ) const
   {
      vertexDoFFunction_.interpolate( expr, level, flag );
      edgeDoFFunction_.interpolate( expr, level, flag );
   }

   inline void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, BoundaryUID boundaryUID ) const
   {
      vertexDoFFunction_.interpolate( expr, level, boundaryUID );
      edgeDoFFunction_.interpolate( expr, level, boundaryUID );
   }

   inline void interpolateExtended( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                                    const std::vector< P2Function< ValueType >* > srcFunctions,
                                    uint_t                                        level,
                                    DoFType                                       flag = All ) const
   {
      std::vector< vertexdof::VertexDoFFunction< ValueType >* > vertexDoFFunctions;
      std::vector< EdgeDoFFunction< ValueType >* >              edgeDoFFunctions;

      for( const auto& function : srcFunctions )
      {
         vertexDoFFunctions.push_back( function->vertexDoFFunction_.get() );
         edgeDoFFunctions.push_back( function->edgeDoFFunction_.get() );
      }

      vertexDoFFunction_.interpolateExtended( expr, vertexDoFFunctions, level, flag );
      edgeDoFFunction_.interpolateExtended( expr, edgeDoFFunctions, level, flag );
   }

   inline void swap( const P2Function< ValueType > & other,
                     const uint_t & level,
                     const DoFType & dofType = All) const
   {
      vertexDoFFunction_.swap( other.getVertexDoFFunction(), level, dofType );
      edgeDoFFunction_.swap( other.getEdgeDoFFunction(), level, dofType );
   }

   inline void assign( const std::vector< ValueType >&                                               scalars,
                       const std::vector< std::reference_wrapper< const P2Function< ValueType > > >& functions,
                       uint_t                                                                        level,
                       DoFType                                                                       flag = All ) const
   {
      std::vector< std::reference_wrapper< const vertexdof::VertexDoFFunction< ValueType > > > vertexDoFFunctions;
      std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > > edgeDoFFunctions;

      for( const P2Function< ValueType >& function : functions )
      {
         vertexDoFFunctions.push_back( function.vertexDoFFunction_ );
         edgeDoFFunctions.push_back( function.edgeDoFFunction_ );
      }

      vertexDoFFunction_.assign( scalars, vertexDoFFunctions, level, flag );
      edgeDoFFunction_.assign( scalars, edgeDoFFunctions, level, flag );
   }

   inline void add( const std::vector< ValueType >&                                               scalars,
                    const std::vector< std::reference_wrapper< const P2Function< ValueType > > >& functions,
                    uint_t                                                                        level,
                    DoFType                                                                       flag = All ) const
   {
      std::vector< std::reference_wrapper< const vertexdof::VertexDoFFunction< ValueType > > > vertexDoFFunctions;
      std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > > edgeDoFFunctions;

      for( const P2Function< ValueType >& function : functions )
      {
         vertexDoFFunctions.push_back( function.vertexDoFFunction_ );
         edgeDoFFunctions.push_back( function.edgeDoFFunction_ );
      }

      vertexDoFFunction_.add( scalars, vertexDoFFunctions, level, flag );
      edgeDoFFunction_.add( scalars, edgeDoFFunctions, level, flag );
   }

   inline real_t dotGlobal( const P2Function< ValueType >& rhs, const uint_t level, const DoFType flag = All ) const
   {
      real_t sum = dotLocal( rhs, level, flag );
      this->startTiming( "Dot (reduce)" );
      walberla::mpi::allReduceInplace( sum, walberla::mpi::SUM, walberla::mpi::MPIManager::instance()->comm() );
      this->stopTiming( "Dot (reduce)" );
      return sum;
   }

   inline real_t dotLocal(const P2Function< ValueType >& rhs, const uint_t level, const DoFType flag = All ) const
   {
      real_t sum = real_c( 0 );
      sum += vertexDoFFunction_.dotLocal( rhs.vertexDoFFunction_, level, flag );
      sum += edgeDoFFunction_.dotLocal( rhs.edgeDoFFunction_, level, flag );
      return sum;
   }

   inline void prolongateP1ToP2( const hhg::P1Function< ValueType >& p1Function,
                                 const uint_t& level,
                                 const DoFType& flag = All ) const
   {
      // Note: 'this' is the dst function - therefore we test this' boundary conditions

      this->startTiming( "Prolongate P1 -> P2" );

      p1Function.template startCommunication< Vertex, Edge >( level );
      p1Function.template startCommunication< Edge, Face >( level );

      for( const auto& it : this->getStorage()->getVertices() )
      {
         const Vertex& vertex = *it.second;

         const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if( testFlag( vertexBC, flag ) )
         {
            P2::macrovertex::prolongateP1ToP2< ValueType >( level,
                                                            vertex,
                                                            vertexDoFFunction_.getVertexDataID(),
                                                            edgeDoFFunction_.getVertexDataID(),
                                                            p1Function.getVertexDataID() );
         }
      }

      p1Function.template endCommunication< Vertex, Edge >( level );

      for( const auto& it : this->getStorage()->getEdges() )
      {
         const Edge& edge = *it.second;

         const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if( testFlag( edgeBC, flag ) )
         {
            P2::macroedge::prolongateP1ToP2< ValueType >( level,
                                                          edge,
                                                          vertexDoFFunction_.getEdgeDataID(),
                                                          edgeDoFFunction_.getEdgeDataID(),
                                                          p1Function.getEdgeDataID() );
         }
      }

      p1Function.template endCommunication< Edge, Face >( level );

      for( const auto& it : this->getStorage()->getFaces() )
      {
         const Face& face = *it.second;

         const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if( testFlag( faceBC, flag ) )
         {
            P2::macroface::prolongateP1ToP2< ValueType >( level,
                                                          face,
                                                          vertexDoFFunction_.getFaceDataID(),
                                                          edgeDoFFunction_.getFaceDataID(),
                                                          p1Function.getFaceDataID() );
         }
      }

      this->stopTiming( "Prolongate P1 -> P2" );
   }

   inline void restrictP2ToP1( const P1Function< ValueType >& p1Function,
                               const uint_t&                                     level,
                               const DoFType&                                    flag = All ) const
   {
      this->startTiming( "Restrict P2 -> P1" );

      vertexDoFFunction_.template startCommunication< Edge, Vertex >( level );
      edgeDoFFunction_.template startCommunication< Edge, Vertex >( level );

      vertexDoFFunction_.template startCommunication< Vertex, Edge >( level );
      edgeDoFFunction_.template startCommunication< Vertex, Edge >( level );

      vertexDoFFunction_.template startCommunication< Face, Edge >( level );
      edgeDoFFunction_.template startCommunication< Face, Edge >( level );

      vertexDoFFunction_.template startCommunication< Edge, Face >( level );
      edgeDoFFunction_.template startCommunication< Edge, Face >( level );

      vertexDoFFunction_.template endCommunication< Edge, Vertex >( level );
      edgeDoFFunction_.template endCommunication< Edge, Vertex >( level );

      for( const auto& it : this->getStorage()->getVertices() )
      {
         const Vertex& vertex = *it.second;

         const DoFType vertexBC = p1Function.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if( testFlag( vertexBC, flag ) )
         {
            P2::macrovertex::restrictP2ToP1< ValueType >( level,
                                                          vertex,
                                                          vertexDoFFunction_.getVertexDataID(),
                                                          edgeDoFFunction_.getVertexDataID(),
                                                          p1Function.getVertexDataID() );
         }
      }

      vertexDoFFunction_.template endCommunication< Vertex, Edge >( level );
      edgeDoFFunction_.template endCommunication< Vertex, Edge >( level );

      vertexDoFFunction_.template endCommunication< Face, Edge >( level );
      edgeDoFFunction_.template endCommunication< Face, Edge >( level );

      for( const auto& it : this->getStorage()->getEdges() )
      {
         const Edge& edge = *it.second;

         const DoFType edgeBC = p1Function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if( testFlag( edgeBC, flag ) )
         {
            P2::macroedge::restrictP2ToP1< ValueType >( level,
                                                        edge,
                                                        vertexDoFFunction_.getEdgeDataID(),
                                                        edgeDoFFunction_.getEdgeDataID(),
                                                        p1Function.getEdgeDataID() );
         }
      }

      vertexDoFFunction_.template endCommunication< Edge, Face >( level );
      edgeDoFFunction_.template endCommunication< Edge, Face >( level );

      for( const auto& it : this->getStorage()->getFaces() )
      {
         const Face& face = *it.second;

         const DoFType faceBC = p1Function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if( testFlag( faceBC, flag ) )
         {
            P2::macroface::restrictP2ToP1< ValueType >( level,
                                                        face,
                                                        vertexDoFFunction_.getFaceDataID(),
                                                        edgeDoFFunction_.getFaceDataID(),
                                                        p1Function.getFaceDataID() );
         }
      }

      this->stopTiming( "Restrict P2 -> P1" );
   }

   inline void restrictInjection( uint_t sourceLevel, DoFType flag = All ) const
   {
      for( const auto& it : this->getStorage()->getFaces() )
      {
         const Face& face = *it.second;

         const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if( testFlag( faceBC, flag ) )
         {
            P2::macroface::restrictInjection< ValueType >(
                sourceLevel, face, vertexDoFFunction_.getFaceDataID(), edgeDoFFunction_.getFaceDataID() );
         }
      }

      for( const auto& it : this->getStorage()->getEdges() )
      {
         const Edge& edge = *it.second;

         const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if( testFlag( edgeBC, flag ) )
         {
            P2::macroedge::restrictInjection< ValueType >(
                sourceLevel, edge, vertexDoFFunction_.getEdgeDataID(), edgeDoFFunction_.getEdgeDataID() );
         }
      }

      for( const auto& it : this->getStorage()->getVertices() )
      {
         const Vertex& vertex = *it.second;

         const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if( testFlag( vertexBC, flag ) )
         {
            P2::macrovertex::restrictInjection< ValueType >(
                sourceLevel, vertex, vertexDoFFunction_.getVertexDataID(), edgeDoFFunction_.getVertexDataID() );
         }
      }
   }

   inline void interpolate( std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                            const std::vector< P2Function< ValueType >* >                                  srcFunctions,
                            uint_t                                                                         level,
                            DoFType                                                                        flag = All ) const
   {
      std::vector< vertexdof::VertexDoFFunction< ValueType >* > vertexDoFFunctions;
      std::vector< EdgeDoFFunction< ValueType >* >              edgeDoFFunctions;

      for( const auto& function : srcFunctions )
      {
         vertexDoFFunctions.push_back( function->vertexDoFFunction_.get() );
         edgeDoFFunctions.push_back( function->edgeDoFFunction_.get() );
      }

      vertexDoFFunction_.interpolateExtended( expr, vertexDoFFunctions, level, flag );
      edgeDoFFunction_.interpolateExtended( expr, edgeDoFFunctions, level, flag );
   }

   inline void prolongateQuadratic( uint_t sourceLevel, DoFType flag = All )
   {
      WALBERLA_ABORT( "P2Function - Prolongate (quadratic) not implemented!" );
   }

   inline void prolongate( uint_t sourceLevel, DoFType flag = All ) const
   {
      edgeDoFFunction_.template communicate< Vertex, Edge >( sourceLevel );
      edgeDoFFunction_.template communicate< Edge, Face >( sourceLevel );

      vertexDoFFunction_.template communicate< Vertex, Edge >( sourceLevel );
      vertexDoFFunction_.template communicate< Edge, Face >( sourceLevel );

      for( const auto& it : this->getStorage()->getFaces() )
      {
         const Face& face = *it.second;

         const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if( testFlag( faceBC, flag ) )
         {
            P2::macroface::prolongate< ValueType >(
                sourceLevel, face, vertexDoFFunction_.getFaceDataID(), edgeDoFFunction_.getFaceDataID() );
         }
      }

      for( const auto& it : this->getStorage()->getEdges() )
      {
         const Edge& edge = *it.second;

         const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if( testFlag( edgeBC, flag ) )
         {
            P2::macroedge::prolongate< ValueType >(
                sourceLevel, edge, vertexDoFFunction_.getEdgeDataID(), edgeDoFFunction_.getEdgeDataID() );
         }
      }

      for( const auto& it : this->getStorage()->getVertices() )
      {
         const Vertex& vertex = *it.second;

         const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if( testFlag( vertexBC, flag ) )
         {
            P2::macrovertex::prolongate< ValueType >(
                sourceLevel, vertex, vertexDoFFunction_.getVertexDataID(), edgeDoFFunction_.getVertexDataID() );
         }
      }
   }

   inline void restrict( uint_t sourceLevel, DoFType flag = All ) const
   {
      edgeDoFFunction_.template communicate< Vertex, Edge >( sourceLevel );
      edgeDoFFunction_.template communicate< Edge, Face >( sourceLevel );

      for( const auto& it : this->getStorage()->getFaces() )
      {
         const Face& face = *it.second;

         const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if( testFlag( faceBC, flag ) )
         {
            P2::macroface::restrict< ValueType >(
                sourceLevel, face, vertexDoFFunction_.getFaceDataID(), edgeDoFFunction_.getFaceDataID() );
         }
      }

      /// sync the vertex dofs which contain the missing edge dofs
      edgeDoFFunction_.template communicate< Face, Edge >( sourceLevel );

      /// remove the temporary updates
      for( const auto& it : this->getStorage()->getFaces() )
      {
         const Face& face = *it.second;

         const DoFType faceBC = this->getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if( testFlag( faceBC, flag ) )
         {
            P2::macroface::postRestrict< ValueType >(
                sourceLevel, face, vertexDoFFunction_.getFaceDataID(), edgeDoFFunction_.getFaceDataID() );
         }
      }

      for( const auto& it : this->getStorage()->getEdges() )
      {
         const Edge& edge = *it.second;

         const DoFType edgeBC = this->getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if( testFlag( edgeBC, flag ) )
         {
            P2::macroedge::restrict< ValueType >(
                sourceLevel, edge, vertexDoFFunction_.getEdgeDataID(), edgeDoFFunction_.getEdgeDataID() );
         }
      }

      //TODO: add real vertex restrict
      for( const auto& it : this->getStorage()->getVertices() )
      {
         const Vertex& vertex = *it.second;

         const DoFType vertexBC = this->getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if( testFlag( vertexBC, flag ) )
         {
            P2::macrovertex::restrictInjection< ValueType >(
                sourceLevel, vertex, vertexDoFFunction_.getVertexDataID(), edgeDoFFunction_.getVertexDataID() );
         }
      }
   }

   inline real_t getMaxValue( uint_t level, DoFType flag = All ) const
   {
      auto localMax = -std::numeric_limits< ValueType >::max();
      localMax      = std::max( localMax, vertexDoFFunction_.getMaxValue( level, flag, false ) );
      localMax      = std::max( localMax, edgeDoFFunction_.getMaxValue( level, flag, false ) );
      walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );

      return localMax;
   }

   inline real_t getMaxMagnitude( uint_t level, DoFType flag = All ) const
   {
      auto localMax = real_t( 0.0 );
      localMax      = std::max( localMax, vertexDoFFunction_.getMaxMagnitude( level, flag, false ) );
      localMax      = std::max( localMax, edgeDoFFunction_.getMaxMagnitude( level, flag, false ) );

      walberla::mpi::allReduceInplace( localMax, walberla::mpi::MAX, walberla::mpi::MPIManager::instance()->comm() );

      return localMax;
   }

   inline real_t getMinValue( uint_t level, DoFType flag = All ) const
   {
      auto localMin = std::numeric_limits< ValueType >::max();
      localMin      = std::min( localMin, vertexDoFFunction_.getMinValue( level, flag, false ) );
      localMin      = std::min( localMin, edgeDoFFunction_.getMinValue( level, flag, false ) );
      walberla::mpi::allReduceInplace( localMin, walberla::mpi::MIN, walberla::mpi::MPIManager::instance()->comm() );

      return localMin;
   }

   inline BoundaryCondition getBoundaryCondition() const
   {
      WALBERLA_ASSERT_EQUAL( vertexDoFFunction_.getBoundaryCondition(),
                             edgeDoFFunction_.getBoundaryCondition(),
                             "P2Function: boundary conditions of underlying vertex- and edgedof functions differ!" );
      return vertexDoFFunction_.getBoundaryCondition();
   }

   inline void enumerate( uint_t level ) const
   {
      this->startTiming( "Enumerate" );

      uint_t counterVertexDoFs = hhg::numberOfLocalDoFs< VertexDoFFunctionTag >( *( this->getStorage() ), level );
      uint_t counterEdgeDoFs   = hhg::numberOfLocalDoFs< EdgeDoFFunctionTag >( *( this->getStorage() ), level );

      std::vector< uint_t > vertexDoFsPerRank = walberla::mpi::allGather( counterVertexDoFs );
      std::vector< uint_t > edgeDoFsPerRank   = walberla::mpi::allGather( counterEdgeDoFs );

      ValueType offset = 0;

      for( uint_t i = 0; i < uint_c( walberla::MPIManager::instance()->rank() ); ++i )
      {
         offset += static_cast< ValueType >( vertexDoFsPerRank[i] + edgeDoFsPerRank[i] );
      }
      enumerate( level, offset );
      this->stopTiming( "Enumerate" );
   }

   inline void enumerate( uint_t level, ValueType& offset ) const
   {
      vertexDoFFunction_.enumerate( level, offset );
      edgeDoFFunction_.enumerate( level, offset );
   }

   inline void setLocalCommunicationMode( const communication::BufferedCommunicator::LocalCommunicationMode& localCommMode )
   {
      vertexDoFFunction_.setLocalCommunicationMode( localCommMode );
      edgeDoFFunction_.setLocalCommunicationMode( localCommMode );
   }

 private:
   using Function< P2Function< ValueType > >::communicators_;

   vertexdof::VertexDoFFunction< ValueType > vertexDoFFunction_;
   EdgeDoFFunction< ValueType >               edgeDoFFunction_;
};

} //namespace hhg
