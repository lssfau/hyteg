#pragma once

#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/p1functionspace/P1Operator.hpp"

namespace hhg
{

class P2P1TaylorHoodStokesOperator
{
public:

  P2P1TaylorHoodStokesOperator(const std::shared_ptr< PrimitiveStorage > & storage, size_t minLevel, size_t maxLevel)
    : A(storage, minLevel, maxLevel),
      div_x_vertexToVertex(storage, minLevel, maxLevel),
      div_x_edgeToVertex(storage, minLevel, maxLevel),
      div_y_vertexToVertex(storage, minLevel, maxLevel),
      div_y_edgeToVertex(storage, minLevel, maxLevel),
      divT_x_vertexToVertex(storage, minLevel, maxLevel),
      divT_x_vertexToEdge(storage, minLevel, maxLevel),
      divT_y_vertexToVertex(storage, minLevel, maxLevel),
      divT_y_vertexToEdge(storage, minLevel, maxLevel)
  {
  }

  void apply(P2P1TaylorHoodFunction<real_t>& src, P2P1TaylorHoodFunction<real_t>& dst, size_t level, DoFType flag)
  {
    A.apply(src.u, dst.u, level, flag, Replace);
    divT_x_vertexToVertex.apply(src.p, *dst.u.getVertexDoFFunction(), level, flag, Add);
    divT_x_vertexToEdge.apply(src.p, *dst.u.getEdgeDoFFunction(), level, flag, Add);

    A.apply(src.v, dst.v, level, flag, Replace);
    divT_y_vertexToVertex.apply(src.p, *dst.v.getVertexDoFFunction(), level, flag, Add);
    divT_y_vertexToEdge.apply(src.p, *dst.v.getEdgeDoFFunction(), level, flag, Add);


    div_x_vertexToVertex.apply(*src.u.getVertexDoFFunction(), dst.p, level, flag | DirichletBoundary, Replace);
    div_x_edgeToVertex.apply(*src.u.getEdgeDoFFunction(), dst.p, level, flag | DirichletBoundary, Add);
    div_y_vertexToVertex.apply(*src.v.getVertexDoFFunction(), dst.p, level, flag | DirichletBoundary, Add);
    div_y_edgeToVertex.apply(*src.v.getEdgeDoFFunction(), dst.p, level, flag | DirichletBoundary, Add);
  }

  P2ConstantLaplaceOperator A;
  P1DivxOperator div_x_vertexToVertex;
  EdgeToVertexDivxOperator div_x_edgeToVertex;
  P1DivyOperator div_y_vertexToVertex;
  EdgeToVertexDivyOperator div_y_edgeToVertex;
  P1DivTxOperator divT_x_vertexToVertex;
  VertexToEdgeDivTxOperator divT_x_vertexToEdge;
  P1DivTyOperator divT_y_vertexToVertex;
  VertexToEdgeDivTyOperator divT_y_vertexToEdge;
};

}
