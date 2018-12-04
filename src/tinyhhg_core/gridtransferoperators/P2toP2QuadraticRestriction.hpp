
#pragma once

#include "tinyhhg_core/p2functionspace/P2Function.hpp"

namespace hhg {

class P2toP2QuadraticRestriction
{
public:

    inline void restrict( const P2Function< real_t > & function, const uint_t & sourceLevel, const DoFType & flag ) const
    {
        const auto storage = function.getStorage();
        const auto& vertexDoFFunction = function.getVertexDoFFunction();
        const auto& edgeDoFFunction = function.getEdgeDoFFunction();
        const auto boundaryCondition = function.getBoundaryCondition();

        edgeDoFFunction.template communicate< Vertex, Edge >( sourceLevel );
        edgeDoFFunction.template communicate< Edge, Face >( sourceLevel );

        for( const auto& it : storage->getFaces() )
        {
            const Face& face = *it.second;

            const DoFType faceBC = boundaryCondition.getBoundaryType( face.getMeshBoundaryFlag() );
            if( testFlag( faceBC, flag ) )
            {
                P2::macroface::restrict< real_t >(
                sourceLevel, face, vertexDoFFunction.getFaceDataID(), edgeDoFFunction.getFaceDataID() );
            }
        }

        /// sync the vertex dofs which contain the missing edge dofs
        edgeDoFFunction.template communicate< Face, Edge >( sourceLevel );

        /// remove the temporary updates
        for( const auto& it : storage->getFaces() )
        {
            const Face& face = *it.second;

            const DoFType faceBC = boundaryCondition.getBoundaryType( face.getMeshBoundaryFlag() );
            if( testFlag( faceBC, flag ) )
            {
                P2::macroface::postRestrict< real_t >(
                sourceLevel, face, vertexDoFFunction.getFaceDataID(), edgeDoFFunction.getFaceDataID() );
            }
        }

        for( const auto& it : storage->getEdges() )
        {
            const Edge& edge = *it.second;

            const DoFType edgeBC = boundaryCondition.getBoundaryType( edge.getMeshBoundaryFlag() );
            if( testFlag( edgeBC, flag ) )
            {
                P2::macroedge::restrict< real_t >(
                sourceLevel, edge, vertexDoFFunction.getEdgeDataID(), edgeDoFFunction.getEdgeDataID() );
            }
        }

        //TODO: add real vertex restrict
        for( const auto& it : storage->getVertices() )
        {
            const Vertex& vertex = *it.second;

            const DoFType vertexBC = boundaryCondition.getBoundaryType( vertex.getMeshBoundaryFlag() );
            if( testFlag( vertexBC, flag ) )
            {
                P2::macrovertex::restrictInjection< real_t >(
                sourceLevel, vertex, vertexDoFFunction.getVertexDataID(), edgeDoFFunction.getVertexDataID() );
            }
        }
    }


};

}