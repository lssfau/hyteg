
#pragma once

#include "tinyhhg_core/p1functionspace/P1Function.hpp"

namespace hhg {

class P1toP1LinearProlongation
{
 public:
   inline void operator()( const P1Function< real_t >& function, const uint_t& sourceLevel, const DoFType& flag ) const
   {
      const uint_t destinationLevel = sourceLevel + 1;

      function.communicate< Vertex, Edge >( sourceLevel );
      function.communicate< Edge, Face >( sourceLevel );
      function.communicate< Face, Edge >( sourceLevel );
      function.communicate< Edge, Vertex >( sourceLevel );

      for( const auto& it : function.getStorage()->getVertices() )
      {
         const Vertex& vertex = *it.second;

         if( testFlag( function.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
         {
            const auto srcData = vertex.getData( function.getVertexDataID() )->getPointer( sourceLevel );
            auto       dstData = vertex.getData( function.getVertexDataID() )->getPointer( destinationLevel );
            prolongateMacroVertex( srcData, dstData, sourceLevel );
         }
      }

      for( const auto& it : function.getStorage()->getEdges() )
      {
         const Edge& edge = *it.second;

         if( testFlag( function.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
         {
            const auto srcData = edge.getData( function.getEdgeDataID() )->getPointer( sourceLevel );
            auto       dstData = edge.getData( function.getEdgeDataID() )->getPointer( destinationLevel );
            prolongateMacroEdge( srcData, dstData, sourceLevel );
         }
      }

      for( const auto& it : function.getStorage()->getFaces() )
      {
         const Face& face = *it.second;

         if( testFlag( function.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
         {
            const auto srcData = face.getData( function.getFaceDataID() )->getPointer( sourceLevel );
            auto       dstData = face.getData( function.getFaceDataID() )->getPointer( destinationLevel );
            prolongateMacroFace( srcData, dstData, sourceLevel );
         }
      }
   }

 private:
   void prolongateMacroVertex( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const;

   void prolongateMacroEdge( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const;

   void prolongateMacroFace( const real_t* src, real_t* dst, const uint_t& sourceLevel ) const;
};

} // namespace hhg