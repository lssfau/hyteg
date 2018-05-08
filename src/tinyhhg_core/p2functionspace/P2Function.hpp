#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/StencilMemory.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFFunction.hpp"
#include "tinyhhg_core/p2functionspace/P2TransferOperators.hpp"

#include "P2Multigrid.hpp"

namespace hhg {

using walberla::real_c;

template < typename ValueType >
class P2Function : public Function< P2Function< ValueType > >
{
 public:
   P2Function( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Function< P2Function< ValueType > >( name, storage, minLevel, maxLevel )
   , vertexDoFFunction_(
         std::make_shared< vertexdof::VertexDoFFunction< ValueType > >( name + "_VertexDoF", storage, minLevel, maxLevel ) )
   , edgeDoFFunction_( std::make_shared< EdgeDoFFunction< ValueType > >( name + "_EdgeDoF", storage, minLevel, maxLevel ) )
   {
      for( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         /// one has to use the communicators of the vertexDoF and edgeDoF function to communicate
         /// TODO: find better solution
         communicators_[level] = NULL;
      }
   }

   std::shared_ptr< vertexdof::VertexDoFFunction< ValueType > > getVertexDoFFunction() const { return vertexDoFFunction_; }
   std::shared_ptr< EdgeDoFFunction< ValueType > >              getEdgeDoFFunction() const { return edgeDoFFunction_; }

   inline void interpolate( std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType flag = All )
   {
      vertexDoFFunction_->interpolate( expr, level, flag );
      edgeDoFFunction_->interpolate( expr, level, flag );
   }

   inline void interpolateExtended( std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                                    const std::vector< P2Function< ValueType >* >                                  srcFunctions,
                                    uint_t                                                                         level,
                                    DoFType                                                                        flag = All )
   {
      std::vector< vertexdof::VertexDoFFunction< ValueType >* > vertexDoFFunctions;
      std::vector< EdgeDoFFunction< ValueType >* >              edgeDoFFunctions;

      for( const auto& function : srcFunctions )
      {
         vertexDoFFunctions.push_back( function->vertexDoFFunction_.get() );
         edgeDoFFunctions.push_back( function->edgeDoFFunction_.get() );
      }

      vertexDoFFunction_->interpolateExtended( expr, vertexDoFFunctions, level, flag );
      edgeDoFFunction_->interpolateExtended( expr, edgeDoFFunctions, level, flag );
   }

   inline void assign( const std::vector< ValueType >                scalars,
                       const std::vector< P2Function< ValueType >* > functions,
                       uint_t                                        level,
                       DoFType                                       flag = All )
   {
      std::vector< vertexdof::VertexDoFFunction< ValueType >* > vertexDoFFunctions;
      std::vector< EdgeDoFFunction< ValueType >* >              edgeDoFFunctions;

      for( const auto& function : functions )
      {
         vertexDoFFunctions.push_back( function->vertexDoFFunction_.get() );
         edgeDoFFunctions.push_back( function->edgeDoFFunction_.get() );
      }

      vertexDoFFunction_->assign( scalars, vertexDoFFunctions, level, flag );
      edgeDoFFunction_->assign( scalars, edgeDoFFunctions, level, flag );
   }

   inline void add( const std::vector< ValueType >                scalars,
                    const std::vector< P2Function< ValueType >* > functions,
                    uint_t                                        level,
                    DoFType                                       flag = All )
   {
      std::vector< vertexdof::VertexDoFFunction< ValueType >* > vertexDoFFunctions;
      std::vector< EdgeDoFFunction< ValueType >* >              edgeDoFFunctions;

      for( const auto& function : functions )
      {
         vertexDoFFunctions.push_back( function->vertexDoFFunction_.get() );
         edgeDoFFunctions.push_back( function->edgeDoFFunction_.get() );
      }

      vertexDoFFunction_->add( scalars, vertexDoFFunctions, level, flag );
      edgeDoFFunction_->add( scalars, edgeDoFFunctions, level, flag );
   }

   inline real_t dot( P2Function< ValueType >& rhs, uint_t level, DoFType flag = All )
   {
      real_t sum = real_c( 0 );
      sum += vertexDoFFunction_->dot( *rhs.vertexDoFFunction_, level, flag );
      sum += edgeDoFFunction_->dot( *rhs.edgeDoFFunction_, level, flag );
      return sum;
   }

   inline void prolongateP1ToP2( const std::shared_ptr< P1Function< ValueType > >& p1Function,
                                 const uint_t&                                     level,
                                 const DoFType&                                    flag = All )
   {
      this->startTiming( "Prolongate P1 -> P2" );

      p1Function->getCommunicator( level )->template startCommunication< Vertex, Edge >();
      p1Function->getCommunicator( level )->template startCommunication< Edge, Face >();

      for( const auto& it : this->getStorage()->getVertices() )
      {
         const Vertex& vertex = *it.second;

         if( testFlag( vertex.getDoFType(), flag ) )
         {
            P2::macrovertex::prolongateP1ToP2< ValueType >( level,
                                                            vertex,
                                                            vertexDoFFunction_->getVertexDataID(),
                                                            edgeDoFFunction_->getVertexDataID(),
                                                            p1Function->getVertexDataID() );
         }
      }

      p1Function->getCommunicator( level )->template endCommunication< Vertex, Edge >();

      for( const auto& it : this->getStorage()->getEdges() )
      {
         const Edge& edge = *it.second;

         if( testFlag( edge.getDoFType(), flag ) )
         {
            P2::macroedge::prolongateP1ToP2< ValueType >( level,
                                                          edge,
                                                          vertexDoFFunction_->getEdgeDataID(),
                                                          edgeDoFFunction_->getEdgeDataID(),
                                                          p1Function->getEdgeDataID() );
         }
      }

      p1Function->getCommunicator( level )->template endCommunication< Edge, Face >();

      for( const auto& it : this->getStorage()->getFaces() )
      {
         const Face& face = *it.second;

         if( testFlag( face.type, flag ) )
         {
            P2::macroface::prolongateP1ToP2< ValueType >( level,
                                                          face,
                                                          vertexDoFFunction_->getFaceDataID(),
                                                          edgeDoFFunction_->getFaceDataID(),
                                                          p1Function->getFaceDataID() );
         }
      }

      this->stopTiming( "Prolongate P1 -> P2" );
   }

   inline void restrictP2ToP1( const std::shared_ptr< P1Function< ValueType > >& p1Function,
                               const uint_t&                                     level,
                               const DoFType&                                    flag = All )
   {
      this->startTiming( "Restrict P2 -> P1" );

      vertexDoFFunction_->getCommunicator( level )->template startCommunication< Edge, Vertex >();
      edgeDoFFunction_->getCommunicator( level )->template startCommunication< Edge, Vertex >();

      vertexDoFFunction_->getCommunicator( level )->template startCommunication< Vertex, Edge >();
      edgeDoFFunction_->getCommunicator( level )->template startCommunication< Vertex, Edge >();

      vertexDoFFunction_->getCommunicator( level )->template startCommunication< Face, Edge >();
      edgeDoFFunction_->getCommunicator( level )->template startCommunication< Face, Edge >();

      vertexDoFFunction_->getCommunicator( level )->template startCommunication< Edge, Face >();
      edgeDoFFunction_->getCommunicator( level )->template startCommunication< Edge, Face >();

      vertexDoFFunction_->getCommunicator( level )->template endCommunication< Edge, Vertex >();
      edgeDoFFunction_->getCommunicator( level )->template endCommunication< Edge, Vertex >();

      for( const auto& it : this->getStorage()->getVertices() )
      {
         const Vertex& vertex = *it.second;

         if( testFlag( vertex.getDoFType(), flag ) )
         {
            P2::macrovertex::restrictP2ToP1< ValueType >( level,
                                                          vertex,
                                                          vertexDoFFunction_->getVertexDataID(),
                                                          edgeDoFFunction_->getVertexDataID(),
                                                          p1Function->getVertexDataID() );
         }
      }

      vertexDoFFunction_->getCommunicator( level )->template endCommunication< Vertex, Edge >();
      edgeDoFFunction_->getCommunicator( level )->template endCommunication< Vertex, Edge >();

      vertexDoFFunction_->getCommunicator( level )->template endCommunication< Face, Edge >();
      edgeDoFFunction_->getCommunicator( level )->template endCommunication< Face, Edge >();

      for( const auto& it : this->getStorage()->getEdges() )
      {
         const Edge& edge = *it.second;

         if( testFlag( edge.getDoFType(), flag ) )
         {
            P2::macroedge::restrictP2ToP1< ValueType >( level,
                                                        edge,
                                                        vertexDoFFunction_->getEdgeDataID(),
                                                        edgeDoFFunction_->getEdgeDataID(),
                                                        p1Function->getEdgeDataID() );
         }
      }

      vertexDoFFunction_->getCommunicator( level )->template endCommunication< Edge, Face >();
      edgeDoFFunction_->getCommunicator( level )->template endCommunication< Edge, Face >();

      for( const auto& it : this->getStorage()->getFaces() )
      {
         const Face& face = *it.second;

         if( testFlag( face.type, flag ) )
         {
            P2::macroface::restrictP2ToP1< ValueType >( level,
                                                        face,
                                                        vertexDoFFunction_->getFaceDataID(),
                                                        edgeDoFFunction_->getFaceDataID(),
                                                        p1Function->getFaceDataID() );
         }
      }

      this->stopTiming( "Restrict P2 -> P1" );
   }

   inline void restrictInjection( uint_t sourceLevel, DoFType flag = All )
   {
      for( const auto& it : this->getStorage()->getFaces() )
      {
         const Face& face = *it.second;

         if( testFlag( face.type, flag ) )
         {
            P2::macroface::restrictInjection< ValueType >(
                sourceLevel, face, vertexDoFFunction_->getFaceDataID(), edgeDoFFunction_->getFaceDataID() );
         }
      }

      for( const auto& it : this->getStorage()->getEdges() )
      {
         const Edge& edge = *it.second;

         if( testFlag( edge.getDoFType(), flag ) )
         {
            P2::macroedge::restrictInjection< ValueType >(
                sourceLevel, edge, vertexDoFFunction_->getEdgeDataID(), edgeDoFFunction_->getEdgeDataID() );
         }
      }

      for( const auto& it : this->getStorage()->getVertices() )
      {
         const Vertex& vertex = *it.second;

         if( testFlag( vertex.getDoFType(), flag ) )
         {
            P2::macrovertex::restrictInjection< ValueType >(
                sourceLevel, vertex, vertexDoFFunction_->getVertexDataID(), edgeDoFFunction_->getVertexDataID() );
         }
      }
   }

   inline void interpolate( std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                            const std::vector< P2Function< ValueType >* >                                  srcFunctions,
                            uint_t                                                                         level,
                            DoFType                                                                        flag = All )
   {
      std::vector< vertexdof::VertexDoFFunction< ValueType >* > vertexDoFFunctions;
      std::vector< EdgeDoFFunction< ValueType >* >              edgeDoFFunctions;

      for( const auto& function : srcFunctions )
      {
         vertexDoFFunctions.push_back( function->vertexDoFFunction_.get() );
         edgeDoFFunctions.push_back( function->edgeDoFFunction_.get() );
      }

      vertexDoFFunction_->interpolateExtended( expr, vertexDoFFunctions, level, flag );
      edgeDoFFunction_->interpolateExtended( expr, edgeDoFFunctions, level, flag );
   }

   inline void prolongateQuadratic( uint_t sourceLevel, DoFType flag = All )
   {
      WALBERLA_ABORT( "P2Function - Prolongate (quadratic) not implemented!" );
   }

   inline void prolongate( uint_t sourceLevel, DoFType flag = All )
   {
      edgeDoFFunction_->getCommunicator( sourceLevel )->template communicate< Vertex, Edge >();
      edgeDoFFunction_->getCommunicator( sourceLevel )->template communicate< Edge, Face >();

      vertexDoFFunction_->getCommunicator( sourceLevel )->template communicate< Vertex, Edge >();
      vertexDoFFunction_->getCommunicator( sourceLevel )->template communicate< Edge, Face >();

      for( const auto& it : this->getStorage()->getFaces() )
      {
         const Face& face = *it.second;

         if( testFlag( face.type, flag ) )
         {
            P2::macroface::prolongate< ValueType >(
                sourceLevel, face, vertexDoFFunction_->getFaceDataID(), edgeDoFFunction_->getFaceDataID() );
         }
      }

      for( const auto& it : this->getStorage()->getEdges() )
      {
         const Edge& edge = *it.second;

         if( testFlag( edge.getDoFType(), flag ) )
         {
            P2::macroedge::prolongate< ValueType >(
                sourceLevel, edge, vertexDoFFunction_->getEdgeDataID(), edgeDoFFunction_->getEdgeDataID() );
         }
      }

      for( const auto& it : this->getStorage()->getVertices() )
      {
         const Vertex& vertex = *it.second;

         if( testFlag( vertex.getDoFType(), flag ) )
         {
            P2::macrovertex::prolongate< ValueType >(
                sourceLevel, vertex, vertexDoFFunction_->getVertexDataID(), edgeDoFFunction_->getVertexDataID() );
         }
      }
   }

   inline void restrict( uint_t sourceLevel, DoFType flag = All )
   {
      edgeDoFFunction_->getCommunicator( sourceLevel )->template communicate< Vertex, Edge >();
      edgeDoFFunction_->getCommunicator( sourceLevel )->template communicate< Edge, Face >();

      for( const auto& it : this->getStorage()->getFaces() )
      {
         const Face& face = *it.second;

         if( testFlag( face.type, flag ) )
         {
            P2::macroface::restrict< ValueType >(
                sourceLevel, face, vertexDoFFunction_->getFaceDataID(), edgeDoFFunction_->getFaceDataID() );
         }
      }

      /// sync the vertex dofs which contain the missing edge dofs
      edgeDoFFunction_->getCommunicator( sourceLevel )->template communicate< Face, Edge >();

      /// remove the temporary updates
      for( const auto& it : this->getStorage()->getFaces() )
      {
         const Face& face = *it.second;

         if( testFlag( face.type, flag ) )
         {
            P2::macroface::postRestrict< ValueType >(
                sourceLevel, face, vertexDoFFunction_->getFaceDataID(), edgeDoFFunction_->getFaceDataID() );
         }
      }

      for( const auto& it : this->getStorage()->getEdges() )
      {
         const Edge& edge = *it.second;

         if( testFlag( edge.getDoFType(), flag ) )
         {
            P2::macroedge::restrict< ValueType >(
                sourceLevel, edge, vertexDoFFunction_->getEdgeDataID(), edgeDoFFunction_->getEdgeDataID() );
         }
      }

      //TODO: add real vertex restrict
      for( const auto& it : this->getStorage()->getVertices() )
      {
         const Vertex& vertex = *it.second;

         if( testFlag( vertex.getDoFType(), flag ) )
         {
            P2::macrovertex::restrictInjection< ValueType >(
                sourceLevel, vertex, vertexDoFFunction_->getVertexDataID(), edgeDoFFunction_->getVertexDataID() );
         }
      }
   }

  inline real_t getMaxMagnitude( uint_t level, DoFType flag = All )
  {
    real_t localMax = real_t(0.0);
    localMax = std::max( localMax, vertexDoFFunction_->getMaxMagnitude( level, flag ) );
    localMax = std::max( localMax, edgeDoFFunction_->getMaxMagnitude( level, flag ) );
    return localMax;
  }


 private:
   using Function< P2Function< ValueType > >::communicators_;
   inline void enumerate_impl( uint_t level, uint_t& num )
   {
      vertexDoFFunction_->enumerate( level, num );
      edgeDoFFunction_->enumerate( level, num );
   }

   std::shared_ptr< vertexdof::VertexDoFFunction< ValueType > > vertexDoFFunction_;
   std::shared_ptr< EdgeDoFFunction< ValueType > >              edgeDoFFunction_;
};

} //namespace hhg
