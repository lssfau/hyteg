
#include "tinyhhg_core/gridtransferoperators/P2toP2QuadraticRestriction.hpp"

#include "tinyhhg_core/gridtransferoperators/generatedKernels/GeneratedKernelsP2MacroFace2D.hpp"

namespace hhg {

using indexing::Index;
using indexing::IndexIncrement;

void P2toP2QuadraticRestriction::restrictAdditively( const P2Function< real_t >& function,
                                                     const uint_t&               sourceLevel,
                                                     const DoFType& ) const
{
   const auto storage = function.getStorage();

   const uint_t fineLevel   = sourceLevel;
   const uint_t coarseLevel = sourceLevel - 1;

   function.communicate< Vertex, Edge >( fineLevel );
   function.communicate< Edge, Face >( fineLevel );

   for ( const auto& faceIt : function.getStorage()->getFaces() )
   {
      const auto face = faceIt.second;

      const auto vertexFineData = face->getData( function.getVertexDoFFunction().getFaceDataID() )->getPointer( fineLevel );
      const auto edgeFineData   = face->getData( function.getEdgeDoFFunction().getFaceDataID() )->getPointer( fineLevel );

      auto vertexCoarseData = face->getData( function.getVertexDoFFunction().getFaceDataID() )->getPointer( coarseLevel );
      auto edgeCoarseData   = face->getData( function.getEdgeDoFFunction().getFaceDataID() )->getPointer( coarseLevel );

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

      typedef edgedof::EdgeDoFOrientation eo;
      std::map< eo, uint_t >              firstIdxFine;
      std::map< eo, uint_t >              firstIdxCoarse;
      for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
      {
        firstIdxFine[e]   = edgedof::macroface::index( fineLevel, 0, 0, e );
        firstIdxCoarse[e] = edgedof::macroface::index( coarseLevel, 0, 0, e );
      }

      P2::macroface::generated::restrict_2D_macroface_P2_update_vertexdofs( &edgeFineData[firstIdxFine[eo::X]],
                                                                            &edgeFineData[firstIdxFine[eo::XY]],
                                                                            &edgeFineData[firstIdxFine[eo::Y]],
                                                                            vertexCoarseData,
                                                                            vertexFineData,
                                                                            static_cast< int64_t >( coarseLevel ),
                                                                            numNeighborFacesEdge0,
                                                                            numNeighborFacesEdge1,
                                                                            numNeighborFacesEdge2,
                                                                            numNeighborFacesVertex0,
                                                                            numNeighborFacesVertex1,
                                                                            numNeighborFacesVertex2 );

      P2::macroface::generated::restrict_2D_macroface_P2_update_edgedofs( &edgeCoarseData[firstIdxCoarse[eo::X]],
                                                                          &edgeCoarseData[firstIdxCoarse[eo::XY]],
                                                                          &edgeCoarseData[firstIdxCoarse[eo::Y]],
                                                                          &edgeFineData[firstIdxFine[eo::X]],
                                                                          &edgeFineData[firstIdxFine[eo::XY]],
                                                                          &edgeFineData[firstIdxFine[eo::Y]],
                                                                          vertexFineData,
                                                                          static_cast< int64_t >( coarseLevel ),
                                                                          numNeighborFacesEdge0,
                                                                          numNeighborFacesEdge1,
                                                                          numNeighborFacesEdge2 );
   }

   function.getVertexDoFFunction().communicateAdditively< Face, Edge >( coarseLevel );
   function.getVertexDoFFunction().communicateAdditively< Face, Vertex >( coarseLevel );

   function.getEdgeDoFFunction().communicateAdditively< Face, Edge >( coarseLevel );
}

void P2toP2QuadraticRestriction::restrictWithPostCommunication( const hhg::P2Function< walberla::real_t >& function,
                                                                const uint_t&                              sourceLevel,
                                                                const hhg::DoFType&                        flag ) const
{
   const auto  storage           = function.getStorage();
   const auto& vertexDoFFunction = function.getVertexDoFFunction();
   const auto& edgeDoFFunction   = function.getEdgeDoFFunction();
   const auto  boundaryCondition = function.getBoundaryCondition();

   edgeDoFFunction.template communicate< Vertex, Edge >( sourceLevel );
   edgeDoFFunction.template communicate< Edge, Face >( sourceLevel );

   for ( const auto& it : storage->getFaces() )
   {
      const Face& face = *it.second;

      const DoFType faceBC = boundaryCondition.getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         P2::macroface::restrict< real_t >(
             sourceLevel, face, vertexDoFFunction.getFaceDataID(), edgeDoFFunction.getFaceDataID() );
      }
   }

   /// sync the vertex dofs which contain the missing edge dofs
   edgeDoFFunction.template communicate< Face, Edge >( sourceLevel );

   /// remove the temporary updates
   for ( const auto& it : storage->getFaces() )
   {
      const Face& face = *it.second;

      const DoFType faceBC = boundaryCondition.getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         P2::macroface::postRestrict< real_t >(
             sourceLevel, face, vertexDoFFunction.getFaceDataID(), edgeDoFFunction.getFaceDataID() );
      }
   }

   for ( const auto& it : storage->getEdges() )
   {
      const Edge& edge = *it.second;

      const DoFType edgeBC = boundaryCondition.getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         P2::macroedge::restrict< real_t >(
             sourceLevel, edge, vertexDoFFunction.getEdgeDataID(), edgeDoFFunction.getEdgeDataID() );
      }
   }

   //TODO: add real vertex restrict
   for ( const auto& it : storage->getVertices() )
   {
      const Vertex& vertex = *it.second;

      const DoFType vertexBC = boundaryCondition.getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         P2::macrovertex::restrictInjection< real_t >(
             sourceLevel, vertex, vertexDoFFunction.getVertexDataID(), edgeDoFFunction.getVertexDataID() );
      }
   }
}

} // namespace hhg
