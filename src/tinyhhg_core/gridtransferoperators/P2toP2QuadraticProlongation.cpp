
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticProlongation.hpp"

#include "tinyhhg_core/gridtransferoperators/generatedKernels/GeneratedKernelsP2MacroFace2D.hpp"

namespace hhg {

void P2toP2QuadraticProlongation::prolongateAdditively( const P2Function< real_t >& function,
                                                        const uint_t&               sourceLevel,
                                                        const DoFType&              flag ) const
{
   WALBERLA_CHECK_EQUAL( function.getVertexDoFFunction().getBoundaryTypeToSkipDuringAdditiveCommunication(),
                         function.getEdgeDoFFunction().getBoundaryTypeToSkipDuringAdditiveCommunication() );
   WALBERLA_CHECK_EQUAL( flag, function.getVertexDoFFunction().getBoundaryTypeToSkipDuringAdditiveCommunication() ^ All );
   WALBERLA_CHECK_EQUAL( flag, function.getEdgeDoFFunction().getBoundaryTypeToSkipDuringAdditiveCommunication() ^ All );

   const auto storage = function.getStorage();

   const uint_t fineLevel   = sourceLevel + 1;
   const uint_t coarseLevel = sourceLevel;

   function.communicate< Vertex, Edge >( coarseLevel );
   function.communicate< Edge, Face >( coarseLevel );

   for ( const auto& faceIt : function.getStorage()->getFaces() )
   {
      const auto face = faceIt.second;

      auto vertexFineData = face->getData( function.getVertexDoFFunction().getFaceDataID() )->getPointer( fineLevel );
      auto edgeFineData   = face->getData( function.getEdgeDoFFunction().getFaceDataID() )->getPointer( fineLevel );

      const auto vertexCoarseData = face->getData( function.getVertexDoFFunction().getFaceDataID() )->getPointer( coarseLevel );
      const auto edgeCoarseData   = face->getData( function.getEdgeDoFFunction().getFaceDataID() )->getPointer( coarseLevel );

      // we need to set the face ghost-layers to zero explicitly since this is not necessarily done by interpolation
      for ( const auto & it : vertexdof::macroface::Iterator( fineLevel, 0 ) )
      {
         vertexFineData[ vertexdof::macroface::index( fineLevel, it.x(), it.y() ) ] = real_c( 0 );
      }
      for ( const auto & it : edgedof::macroface::Iterator( fineLevel, 0 ) )
      {
         edgeFineData[ edgedof::macroface::index( fineLevel, it.x(), it.y(), edgedof::EdgeDoFOrientation::X ) ] = real_c( 0 );
         edgeFineData[ edgedof::macroface::index( fineLevel, it.x(), it.y(), edgedof::EdgeDoFOrientation::XY ) ] = real_c( 0 );
         edgeFineData[ edgedof::macroface::index( fineLevel, it.x(), it.y(), edgedof::EdgeDoFOrientation::Y ) ] = real_c( 0 );
      }

      const double numNeighborFacesEdge0 =
          static_cast< double >( storage->getEdge( face->neighborEdges().at( 0 ) )->getNumNeighborFaces() );
      const double numNeighborFacesEdge1 =
          static_cast< double >( storage->getEdge( face->neighborEdges().at( 1 ) )->getNumNeighborFaces() );
      const double numNeighborFacesEdge2 =
          static_cast< double >( storage->getEdge( face->neighborEdges().at( 2 ) )->getNumNeighborFaces() );
      const double numNeighborFacesVertex0 =
          static_cast< double >( storage->getVertex( face->neighborVertices().at( 0 ) )->getNumNeighborFaces() );
      const double numNeighborFacesVertex1 =
          static_cast< double >( storage->getVertex( face->neighborVertices().at( 1 ) )->getNumNeighborFaces() );
      const double numNeighborFacesVertex2 =
          static_cast< double >( storage->getVertex( face->neighborVertices().at( 2 ) )->getNumNeighborFaces() );

      P2::macroface::generated::prolongate_2D_macroface_P2_push_from_vertexdofs( edgeFineData,
                                                                                 vertexCoarseData,
                                                                                 vertexFineData,
                                                                                 static_cast< int64_t >( coarseLevel ),
                                                                                 numNeighborFacesEdge0,
                                                                                 numNeighborFacesEdge1,
                                                                                 numNeighborFacesEdge2,
                                                                                 numNeighborFacesVertex0,
                                                                                 numNeighborFacesVertex1,
                                                                                 numNeighborFacesVertex2 );

      P2::macroface::generated::prolongate_2D_macroface_P2_push_from_edgedofs( edgeCoarseData,
                                                                               edgeFineData,
                                                                               vertexFineData,
                                                                               static_cast< int64_t >( coarseLevel ),
                                                                               numNeighborFacesEdge0,
                                                                               numNeighborFacesEdge1,
                                                                               numNeighborFacesEdge2 );
   }

   function.getVertexDoFFunction().communicateAdditively< Face, Edge >( fineLevel );
   function.getVertexDoFFunction().communicateAdditively< Face, Vertex >( fineLevel );

   function.getEdgeDoFFunction().communicateAdditively< Face, Edge >( fineLevel );
}

void P2toP2QuadraticProlongation::prolongateStandard( const P2Function< real_t >& function,
                                                      const uint_t&               sourceLevel,
                                                      const DoFType&              flag ) const
{
   const auto& vertexDoFFunction = function.getVertexDoFFunction();
   const auto& edgeDoFFunction   = function.getEdgeDoFFunction();

   edgeDoFFunction.template communicate< Vertex, Edge >( sourceLevel );
   edgeDoFFunction.template communicate< Edge, Face >( sourceLevel );

   vertexDoFFunction.template communicate< Vertex, Edge >( sourceLevel );
   vertexDoFFunction.template communicate< Edge, Face >( sourceLevel );

   for ( const auto& it : function.getStorage()->getFaces() )
   {
      const Face& face = *it.second;

      const DoFType faceBC = function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         P2::macroface::prolongate< real_t >(
             sourceLevel, face, vertexDoFFunction.getFaceDataID(), edgeDoFFunction.getFaceDataID() );
      }
   }

   for ( const auto& it : function.getStorage()->getEdges() )
   {
      const Edge& edge = *it.second;

      const DoFType edgeBC = function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         P2::macroedge::prolongate< real_t >(
             sourceLevel, edge, vertexDoFFunction.getEdgeDataID(), edgeDoFFunction.getEdgeDataID() );
      }
   }

   for ( const auto& it : function.getStorage()->getVertices() )
   {
      const Vertex& vertex = *it.second;

      const DoFType vertexBC = function.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         P2::macrovertex::prolongate< real_t >(
             sourceLevel, vertex, vertexDoFFunction.getVertexDataID(), edgeDoFFunction.getVertexDataID() );
      }
   }
}
} // namespace hhg