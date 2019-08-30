
#pragma once

#include "hyteg/Operator.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"

namespace hyteg {

class P1PointwiseOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
 public:
   P1PointwiseOperator( const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   {
      auto faceP1StencilMemoryDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Face > >(
          minLevel_, maxLevel_, vertexDoFMacroFacePointwiseStencilMemorySize );
      auto edgeP1StencilMemoryDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Edge > >(
          minLevel_, maxLevel_, vertexDoFMacroEdgePointwiseStencilMemorySize );
      auto vertexP1StencilMemoryDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Vertex > >(
          minLevel_, maxLevel_, vertexDoFMacroVertexPointwiseStencilMemorySize );

      storage->addFaceData( faceStencilID_, faceP1StencilMemoryDataHandling, "P1PointwiseOperatorFaceStencil" );
      storage->addEdgeData( edgeStencilID_, edgeP1StencilMemoryDataHandling, "P1PointwiseOperatorEdgeStencil" );
      storage->addVertexData( vertexStencilID_, vertexP1StencilMemoryDataHandling, "P1PointwiseOperatorVertexStencil" );
   }

   ~P1PointwiseOperator() override = default;

   void apply( const P1Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      this->startTiming( "Apply" );
      communication::syncFunctionBetweenPrimitives( src, level );

      for ( const auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if ( testFlag( vertexBC, flag ) )
         {
            vertexdof::macrovertex::applyPointwise< real_t >(
                level, vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), updateType );
         }
      }

      for ( const auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            vertexdof::macroedge::applyPointwise< real_t >(
                level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
         }
      }

      for ( const auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            vertexdof::macroface::applyPointwise< real_t >(
                level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
         }
      }
      this->stopTiming( "Apply" );
   }

   const PrimitiveDataID< StencilMemory< real_t >, Vertex >& getVertexStencilID() const { return vertexStencilID_; }
   const PrimitiveDataID< StencilMemory< real_t >, Edge >&   getEdgeStencilID() const { return edgeStencilID_; }
   const PrimitiveDataID< StencilMemory< real_t >, Face >&   getFaceStencilID() const { return faceStencilID_; }
   const PrimitiveDataID< StencilMemory< real_t >, Cell >&   getCellStencilID() const { return cellStencilID_; }

 private:
   PrimitiveDataID< StencilMemory< real_t >, Vertex > vertexStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Edge >   edgeStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Face >   faceStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Cell >   cellStencilID_;
};

} // namespace hyteg
