#include "VertexDoFIndexing.hpp"

#include "tinyhhg_core/indexing/MacroEdgeIndexing.hpp"

namespace hhg {
namespace vertexdof {

// ##################
// ### Macro Edge ###
// ##################

namespace macroedge {

uint_t neighborFaceGhostLayerSize( const uint_t& level )
{
   return levelinfo::num_microvertices_per_edge( level ) - 1;
}

uint_t neighborCellGhostLayerSize( const uint_t& level )
{
   return levelinfo::num_microvertices_per_edge( level ) - 2;
}

uint_t index( const uint_t& level, const uint_t& x )
{
   return hhg::indexing::macroEdgeIndex( levelinfo::num_microvertices_per_edge( level ), x );
}

uint_t indexOnNeighborFace( const uint_t& level, const uint_t& x, const uint_t& neighbor )
{
   return hhg::indexing::macroEdgeSize( levelinfo::num_microvertices_per_edge( level ) ) +
          neighbor * hhg::indexing::macroEdgeSize( neighborFaceGhostLayerSize( level ) ) +
          hhg::indexing::macroEdgeIndex( neighborFaceGhostLayerSize( level ), x );
}

uint_t indexOnNeighborCell( const uint_t& level,
                                                                 const uint_t& x,
                                                                 const uint_t& neighbor,
                                                                 const uint_t& numNeighborFaces )
{
   return hhg::indexing::macroEdgeSize( levelinfo::num_microvertices_per_edge( level ) ) +
          numNeighborFaces * hhg::indexing::macroEdgeSize( neighborFaceGhostLayerSize( level ) ) +
          neighbor * hhg::indexing::macroEdgeSize( neighborCellGhostLayerSize( level ) ) +
          hhg::indexing::macroEdgeIndex( neighborCellGhostLayerSize( level ), x );
}

uint_t indexFromVertexOnNeighborFace(const uint_t &level, const uint_t &x, const uint_t &faceID, const stencilDirection &dir) {
   typedef stencilDirection sD;
   WALBERLA_ASSERT( dir == sD::VERTEX_W || dir == sD::VERTEX_E );
   switch( dir )
   {
      case sD::VERTEX_W:
         return indexOnNeighborFace( level, x - 1, faceID );
      case sD::VERTEX_E:
         return indexOnNeighborFace( level, x, faceID  );
      default:
         return std::numeric_limits< uint_t >::max();
   }
}

uint_t indexFromVertexOnNeighborCell(const uint_t &level, const uint_t &x, const uint_t &cellID, const uint_t &numNeighborFaces) {
   WALBERLA_ASSERT_GREATER_EQUAL( x, 1, "The 0th edge idx has no cell neighbor." );
   WALBERLA_ASSERT_LESS_EQUAL   ( x, neighborCellGhostLayerSize( level ) );
   return indexOnNeighborCell( level, x - 1, cellID, numNeighborFaces );
}

uint_t indexFromVertex(const uint_t &level, const uint_t &x, const stencilDirection &dir) {
   typedef stencilDirection sD;

   switch( dir )
   {
      case sD::VERTEX_C:
         return index( level, x );
      case sD::VERTEX_E:
         return index( level, x + 1 );
      case sD::VERTEX_W:
         return index( level, x - 1 );
      case sD::VERTEX_N:
         return indexFromVertexOnNeighborFace( level, x, 1, sD::VERTEX_E );
      case sD::VERTEX_S:
         return indexFromVertexOnNeighborFace( level, x, 0, sD::VERTEX_W );
      case sD::VERTEX_NW:
         return indexFromVertexOnNeighborFace( level, x, 1, sD::VERTEX_W );
      case sD::VERTEX_SE:
         return indexFromVertexOnNeighborFace( level, x, 0, sD::VERTEX_E );
      default:
         return std::numeric_limits< uint_t >::max();
   }
}


uint_t stencilIndexOnEdge(const stencilDirection &dir) {
   typedef stencilDirection sD;
   WALBERLA_ASSERT( dir == sD::VERTEX_C || dir == sD::VERTEX_W || dir == sD::VERTEX_E );
   switch ( dir )
   {
      case sD::VERTEX_C:
         return 3;
      case sD::VERTEX_W:
         return 2;
      case sD::VERTEX_E:
         return 4;
      default:
         return std::numeric_limits< uint_t >::max();
   }
}

uint_t stencilIndexOnNeighborFace(const stencilDirection &dir, const uint_t &faceID) {
   typedef stencilDirection sD;
   WALBERLA_ASSERT( dir == sD::VERTEX_W || dir == sD::VERTEX_E );
   if ( faceID == 0 )
   {
      switch ( dir )
      {
         case sD::VERTEX_W:
            return 0;
         case sD::VERTEX_E:
            return 1;
         default:
         WALBERLA_ABORT("wrong direction")
      }
   } else {
      switch (dir) {
         case sD::VERTEX_W:
            return 3 + 2 * faceID;
         case sD::VERTEX_E:
            return 3 + 2 * faceID + 1;
         default:
         WALBERLA_ABORT("wrong direction")
      }
   }
}

uint_t stencilIndexOnNeighborCell(const uint_t &cellID, const uint_t &numNeighborFaces) {
   return 3 + 2 * numNeighborFaces + cellID;
}

uint_t indexFromHorizontalEdge(const uint_t &level, const uint_t &x, const stencilDirection &dir) {
   typedef stencilDirection sD;

   switch( dir )
   {
      case sD::VERTEX_W:
         return index( level, x );
      case sD::VERTEX_E:
         return index( level, x + 1 );
      case sD::VERTEX_SE:
         return indexOnNeighborFace( level, x, 0 );
      case sD::VERTEX_NW:
         return indexOnNeighborFace( level, x, 1 );
      default:
         return std::numeric_limits< uint_t >::max();
   }
}

Iterator::Iterator(const uint_t &level, const uint_t &offsetToCenter) :
   EdgeIterator( levelinfo::num_microvertices_per_edge( level ), offsetToCenter )
{}
} // namespace macroedge

// ##################
// ### Macro Face ###
// ##################

namespace macroface {

uint_t index(const uint_t &level, const uint_t &x, const uint_t &y) {
   return hhg::indexing::macroFaceIndex( levelinfo::num_microvertices_per_edge( level ), x, y );
}

uint_t index(const uint_t &level, const uint_t &x, const uint_t &y, const uint_t &neighbor) {
   assert( neighbor <= 1 );

   return hhg::indexing::macroFaceSize( levelinfo::num_microvertices_per_edge( level ) )
          + neighbor * hhg::indexing::macroFaceSize( levelinfo::num_microvertices_per_edge( level ) - 1 )
          + hhg::indexing::macroFaceIndex( levelinfo::num_microvertices_per_edge( level ) - 1, x, y );

}

uint_t indexFromVertex(const uint_t &level, const uint_t &x, const uint_t &y, const stencilDirection &dir) {
   typedef stencilDirection sD;

   switch( dir )
   {
      case sD::VERTEX_C:
         return index( level, x, y );
      case sD::VERTEX_E:
         return index( level, x + 1, y );
      case sD::VERTEX_W:
         return index( level, x - 1, y );
      case sD::VERTEX_N:
         return index( level, x, y + 1 );
      case sD::VERTEX_S:
         return index( level, x, y - 1 );
      case sD::VERTEX_NW:
         return index( level, x - 1, y + 1 );
      case sD::VERTEX_SE:
         return index( level, x + 1, y - 1 );
      case sD::VERTEX_TC:
         return index( level, x, y, 0 );
      case sD::VERTEX_TW:
         return index( level, x - 1, y, 0 );
      case sD::VERTEX_TS:
         return index( level, x, y - 1, 0 );
      case sD::VERTEX_TSE:
         return index( level, x + 1, y - 1, 0 );
      case sD::VERTEX_TSW:
         return index( level, x - 1, y - 1, 0 );
      case sD::VERTEX_TNW:
         return index( level, x - 1, y + 1, 0 );
      case sD::VERTEX_BC:
         return index( level, x, y, 1 );
      case sD::VERTEX_BN:
         return index( level, x, y + 1, 1 );
      case sD::VERTEX_BS:
         return index( level, x, y - 1, 1 );
      case sD::VERTEX_BE:
         return index( level, x + 1, y, 1 );
      case sD::VERTEX_BW:
         return index( level, x - 1, y, 1 );
      case sD::VERTEX_BSW:
         return index( level, x - 1, y - 1, 1 );
      case sD::VERTEX_BSE:
         return index( level, x + 1, y - 1, 1 );
      case sD::VERTEX_BNW:
         return index( level, x - 1, y + 1, 1 );
      default:
         return std::numeric_limits< uint_t >::max();
   }
}

uint_t indexFromHorizontalEdge(const uint_t &level, const uint_t &x, const uint_t &y, const stencilDirection &dir) {
   typedef stencilDirection sD;

   switch( dir )
   {
      case sD::VERTEX_W:
         return index( level, x, y );
      case sD::VERTEX_E:
         return index( level, x + 1, y );
      case sD::VERTEX_SE:
         return index( level, x + 1, y - 1 );
      case sD::VERTEX_NW:
         return index( level, x, y + 1 );
      default:
         return std::numeric_limits< uint_t >::max();
   }
}

uint_t indexFromDiagonalEdge(const uint_t &level, const uint_t &x, const uint_t &y, const stencilDirection &dir) {
   typedef stencilDirection sD;

   switch( dir )
   {
      case sD::VERTEX_SE:
         return index( level, x + 1, y );
      case sD::VERTEX_NE:
         return index( level, x + 1, y + 1 );
      case sD::VERTEX_NW:
         return index( level, x, y + 1 );
      case sD::VERTEX_SW:
         return index( level, x, y );
      default:
         return std::numeric_limits< uint_t >::max();
   }
}

uint_t indexFromVerticalEdge(const uint_t &level, const uint_t &x, const uint_t &y, const stencilDirection &dir) {
   typedef stencilDirection sD;

   switch( dir )
   {
      case sD::VERTEX_S:
         return index( level, x, y );
      case sD::VERTEX_SE:
         return index( level, x + 1, y );
      case sD::VERTEX_N:
         return index( level, x, y + 1 );
      case sD::VERTEX_NW:
         return index( level, x - 1, y + 1 );
      default:
         return std::numeric_limits< uint_t >::max();
   }
}

uint_t indexFromGrayFace(const uint_t &level, const uint_t &x, const uint_t &y, const stencilDirection &dir) {
   typedef stencilDirection sD;

   switch( dir )
   {
      case sD::VERTEX_SW:
         return index( level, x, y );
      case sD::VERTEX_SE:
         return index( level, x + 1, y );
      case sD::VERTEX_NW:
         return index( level, x, y + 1 );
      default:
         return std::numeric_limits< uint_t >::max();
   }
}

uint_t indexFromBlueFace(const uint_t &level, const uint_t &x, const uint_t &y, const stencilDirection &dir) {
   typedef stencilDirection sD;

   switch( dir )
   {
      case sD::VERTEX_SE:
         return index( level, x + 1, y );
      case sD::VERTEX_NW:
         return index( level, x, y + 1 );
      case sD::VERTEX_NE:
         return index( level, x + 1, y + 1 );
      default:
         return std::numeric_limits< uint_t >::max();
   }
}

bool isVertexOnBoundary(const uint_t &level, const hhg::indexing::Index &idx) {
   if (idx.row() == 0 ){
      return true;
   } else if( idx.col() == 0 ){
      return true;
   } else if( (idx.row() + idx.col()) == ( hhg::levelinfo::num_microvertices_per_edge( level ) - 1 ) ){
      return true;
   } else {
      return false;
   }
}

Iterator::Iterator(const uint_t &level, const uint_t &offsetToCenter) :
   FaceIterator( levelinfo::num_microvertices_per_edge( level ), offsetToCenter )
{}

BorderIterator::BorderIterator(const uint_t &level, const hhg::indexing::FaceBorderDirection &direction, const uint_t &offsetToCenter,
                               const uint_t &offsetFromVertices) :
   FaceBorderIterator( levelinfo::num_microvertices_per_edge( level ), direction, offsetToCenter, offsetFromVertices )
{}
}

} // namespace vertexdof
} // namespace hhg