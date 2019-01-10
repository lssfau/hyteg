#pragma once

#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/p2functionspace/P2Elements.hpp"

#include "tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

#include "generated/p1_to_p2_divt.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

namespace hhg {

using walberla::real_t;

template< class UFCOperator2D, class UFCOperator3D >
class P1ToP2ConstantOperator : public Operator<P1Function < real_t>, P2Function<real_t> > {
 public:

  P1ToP2ConstantOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
      : Operator(storage, minLevel, maxLevel), vertexToVertex(storage, minLevel, maxLevel),
        vertexToEdge(storage, minLevel, maxLevel)
  {
    using namespace P2Elements;

    if ( !storage->hasGlobalCells() )
    {
      // Initialize memory for local 6x6 matrices
      Matrixr< 6, 3 > local_stiffness_gray;
      Matrixr< 6, 3 > local_stiffness_blue;

      // Assemble stencils on all levels
      for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
      {

        // Assemble face stencils
        for ( auto & it : storage_->getFaces())
        {
          Face & face = *it.second;

          // Compute both local stiffness matrices
          compute_local_stiffness( face, level, local_stiffness_gray, fenics::GRAY );
          compute_local_stiffness( face, level, local_stiffness_blue, fenics::BLUE );

//        WALBERLA_LOG_DEVEL_ON_ROOT("local_stiffness_gray =\n" << local_stiffness_gray);
//        WALBERLA_LOG_DEVEL_ON_ROOT("local_stiffness_blue =\n" << local_stiffness_blue);

          // Assemble vertexToVertex stencil
          real_t *vStencil = storage_->getFace( face.getID())->getData( vertexToVertex.getFaceStencilID())->getPointer( level );
          P2Face::VertexToVertex::assembleStencil( local_stiffness_gray, local_stiffness_blue, vStencil );
//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToVertex/Face = {}", PointND<real_t, 7>(&vStencil[0])));

          // Assemble vertexToEdge stencil
          vStencil = storage_->getFace( face.getID())->getData( vertexToEdge.getFaceStencilID())->getPointer( level );
          P2Face::VertexToEdge::assembleStencil( local_stiffness_gray, local_stiffness_blue, vStencil );
//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToEdge/Face = {}", PointND<real_t, 12>(&vStencil[0])));
        }

        // Assemble edge stencils
        for ( auto & it : storage_->getEdges())
        {
          Edge & edge = *it.second;

          // Assemble vertexToVertex stencil
          Face *face = storage_->getFace( edge.neighborFaces()[0] );
          real_t *vStencil = storage_->getEdge( edge.getID())->getData( vertexToVertex.getEdgeStencilID())->getPointer( level );
          compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
          compute_local_stiffness( *face, level, local_stiffness_blue, fenics::BLUE );
          P2Edge::VertexToVertex::assembleStencil( edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, true );

          if ( edge.getNumNeighborFaces() == 2 )
          {
            face = storage_->getFace( edge.neighborFaces()[1] );
            compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
            compute_local_stiffness( *face, level, local_stiffness_blue, fenics::BLUE );
            P2Edge::VertexToVertex::assembleStencil( edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, false );
          }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToVertex/Edge = {}", PointND<real_t, 7>(&vStencil[0])));

          // Assemble vertexToEdge stencil
          face = storage_->getFace( edge.neighborFaces()[0] );
          vStencil = storage_->getEdge( edge.getID())->getData( vertexToEdge.getEdgeStencilID())->getPointer( level );
          compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
          compute_local_stiffness( *face, level, local_stiffness_blue, fenics::BLUE );
          P2Edge::VertexToEdge::assembleStencil( edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, true );

          if ( edge.getNumNeighborFaces() == 2 )
          {
            face = storage_->getFace( edge.neighborFaces()[1] );
            compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
            compute_local_stiffness( *face, level, local_stiffness_blue, fenics::BLUE );
            P2Edge::VertexToEdge::assembleStencil( edge, *face, local_stiffness_gray, local_stiffness_blue, vStencil, false );
          }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToEdge/Edge = {}", PointND<real_t, 4>(&vStencil[0])));
        }

        for ( auto & it : storage_->getVertices())
        {
          Vertex & vertex = *it.second;

          // Assemble VertexToVertex
          real_t *vStencil = storage_->getVertex( vertex.getID())->getData( vertexToVertex.getVertexStencilID())->getPointer( level );
          for ( auto & faceId : vertex.neighborFaces())
          {
            Face *face = storage_->getFace( faceId );
            compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
            P2Vertex::VertexToVertex::assembleStencil( vertex, *face, local_stiffness_gray, vStencil, storage_ );
          }

//        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("vertexToVertex/Vertex = {}", PointND<real_t, 5>(&vStencil[0])));
        }

      }
    }

  }


  P1ConstantOperator< fenics::NoAssemble, UFCOperator3D > const & getVertexToVertexOpr() const {
    return vertexToVertex;
  }

  VertexDoFToEdgeDoFOperator< fenics::NoAssemble, UFCOperator3D > const & getVertexToEdgeOpr() const {
    return vertexToEdge;
  }

   void apply( const P1Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      vertexToVertex.apply( src, dst.getVertexDoFFunction(), level, flag, updateType );
      vertexToEdge.apply( src, dst.getEdgeDoFFunction(), level, flag, updateType );
  }

  void smooth_gs(P1Function< real_t > & dst, P2Function< real_t > & rhs, size_t level, DoFType flag)
  {
    WALBERLA_ABORT("not implemented");
  }

private:

  P1ConstantOperator< fenics::NoAssemble, UFCOperator3D > vertexToVertex;
  VertexDoFToEdgeDoFOperator< fenics::NoAssemble, UFCOperator3D > vertexToEdge;

  void compute_local_stiffness(const Face &face, size_t level, Matrixr<6, 3>& local_stiffness, fenics::ElementType element_type) {
    real_t coords[6];
    fenics::compute_micro_coords(face, level, coords, element_type);
    UFCOperator2D gen;
    gen.tabulate_tensor(local_stiffness.data(), NULL, coords, 0);
  }

};

typedef P1ToP2ConstantOperator< p1_to_p2_divt_cell_integral_0_otherwise, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > P1ToP2ConstantDivTxOperator;
typedef P1ToP2ConstantOperator< p1_to_p2_divt_cell_integral_1_otherwise, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > P1ToP2ConstantDivTyOperator;
typedef P1ToP2ConstantOperator< fenics::NoAssemble,                      p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > P1ToP2ConstantDivTzOperator;

}
