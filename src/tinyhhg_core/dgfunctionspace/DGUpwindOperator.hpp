#pragma once

#include <array>

#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/dgfunctionspace/DGFunction.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

namespace hhg {

template < class VelocityBaseType >
class DGUpwindOperator : public Operator< DGFunction< real_t >, DGFunction< real_t > >
{
   typedef std::array< VelocityBaseType , 2 > VelocityType;

 public:
   DGUpwindOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                     const VelocityType&                        velocity,
                     size_t                                     minLevel,
                     size_t                                     maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , velocity_( velocity )
   {}

   ~DGUpwindOperator() = default;

   void apply( const DGFunction< real_t >& src,
               const DGFunction< real_t >& dst,
               uint_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      // start pulling edge halos
      src.startCommunication< Face, Edge >( level );

      // end pulling edge halos
      src.endCommunication< Face, Edge >( level );

      // start pulling vertex halos
      src.startCommunication< Edge, Vertex >( level );

      // end pulling vertex halos
      src.endCommunication< Edge, Vertex >( level );

      for( auto velocityComponent : velocity_ )
      {
         velocityComponent.template startCommunication< Edge, Vertex >( level );
      }

      for( auto velocityComponent : velocity_ )
      {
         velocityComponent.template startCommunication< Face, Edge >( level );
      }

      for( auto velocityComponent : velocity_ )
      {
         velocityComponent.template endCommunication< Edge, Vertex >( level );
      }

      for( auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if( testFlag( vertexBC, flag ) )
         {
            DGVertex::upwind< real_t >( level,
                                        vertex,
                                        storage_,
                                        src.getVertexDataID(),
                                        dst.getVertexDataID(),
                                        std::array< PrimitiveDataID< FunctionMemory< real_t >, Vertex >, 2 >{
                                            {velocity_[0].getVertexDataID(), velocity_[1].getVertexDataID()}},
                                        updateType );
         }
      }

      dst.startCommunication< Vertex, Edge >( level );

      for( auto velocityComponent : velocity_ )
      {
         velocityComponent.template endCommunication< Face, Edge >( level );
      }

      for( auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if( testFlag( edgeBC, flag ) )
         {
            DGEdge::upwind< real_t >( level,
                                      edge,
                                      storage_,
                                      src.getEdgeDataID(),
                                      dst.getEdgeDataID(),
                                      std::array< PrimitiveDataID< FunctionMemory< real_t >, Edge >, 2 >{
                                          {velocity_[0].getEdgeDataID(), velocity_[1].getEdgeDataID()}},
                                      updateType );
         }
      }

      dst.endCommunication< Vertex, Edge >( level );

      dst.startCommunication< Edge, Face >( level );

      for( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if( testFlag( faceBC, flag ) )
         {
            DGFace::upwind< real_t >( level,
                                      face,
                                      storage_,
                                      src.getFaceDataID(),
                                      dst.getFaceDataID(),
                                      std::array< PrimitiveDataID< FunctionMemory< real_t >, Face >, 2 >{
                                          {velocity_[0].getFaceDataID(), velocity_[1].getFaceDataID()}},
                                      updateType );
         }
      }

      dst.endCommunication< Edge, Face >( level );
   }

 private:
   VelocityType velocity_;
};

} // namespace hhg
