#pragma once

namespace hhg{
namespace stencilWeights{

class tri_1el {
public:
  real_t vertexToVertexStencil[7];
  real_t edgeToVertexStencil[12];
  real_t vertexToEdgeStencil[12];
  real_t edgeToEdgeStencil[15];

  tri_1el() {
    vertexToVertexStencil[0] = 1.0 / 3.0;
    vertexToVertexStencil[1] = 0.0;
    vertexToVertexStencil[2] = 1.0 / 3.0;
    vertexToVertexStencil[3] = 4.0;
    vertexToVertexStencil[4] = 1.0 / 3.0;
    vertexToVertexStencil[5] = 0.0;
    vertexToVertexStencil[6] = 1.0 / 3.0;

    edgeToVertexStencil[0] = -(1.0 + 1.0 / 3.0);
    edgeToVertexStencil[1] = 0.0;
    edgeToVertexStencil[2] = -(1.0 + 1.0 / 3.0);
    edgeToVertexStencil[3] = 0.0;
    edgeToVertexStencil[4] = 0.0;
    edgeToVertexStencil[5] = 0.0;
    edgeToVertexStencil[6] = -(1.0 + 1.0 / 3.0);
    edgeToVertexStencil[7] = 0.0;
    edgeToVertexStencil[8] = -(1.0 + 1.0 / 3.0);
    edgeToVertexStencil[9] = 0.0;
    edgeToVertexStencil[10] = 0.0;
    edgeToVertexStencil[11] = 0.0;

    vertexToEdgeStencil[0] = - (1.0 + 1.0/3.0);
    vertexToEdgeStencil[1] = - (1.0 + 1.0/3.0);
    vertexToEdgeStencil[2] = 0.0;
    vertexToEdgeStencil[3] = 0.0;
    vertexToEdgeStencil[4] = 0.0;
    vertexToEdgeStencil[5] = 0.0;
    vertexToEdgeStencil[6] = 0.0;
    vertexToEdgeStencil[7] = 0.0;
    vertexToEdgeStencil[8] = - (1.0 + 1.0/3.0);
    vertexToEdgeStencil[9] = 0.0;
    vertexToEdgeStencil[10] = - (1.0 + 1.0/3.0);;
    vertexToEdgeStencil[11] = 0.0;

    edgeToEdgeStencil[0] = (5.0 + 1.0/3.0);
    edgeToEdgeStencil[1] = - (1.0 + 1.0/3.0);
    edgeToEdgeStencil[2] = 0.0;
    edgeToEdgeStencil[3] = - (1.0 + 1.0/3.0);
    edgeToEdgeStencil[4] = 0.0;
    edgeToEdgeStencil[5] = (5.0 + 1.0/3.0);
    edgeToEdgeStencil[6] = - (1.0 + 1.0/3.0);
    edgeToEdgeStencil[7] = - (1.0 + 1.0/3.0);
    edgeToEdgeStencil[8] = - (1.0 + 1.0/3.0);
    edgeToEdgeStencil[9] = - (1.0 + 1.0/3.0);
    edgeToEdgeStencil[10] = (5.0 + 1.0/3.0);
    edgeToEdgeStencil[11] = 0.0;
    edgeToEdgeStencil[12] = - (1.0 + 1.0/3.0);
    edgeToEdgeStencil[13] = 0.0;
    edgeToEdgeStencil[14] = - (1.0 + 1.0/3.0);

  }
};


}// namespace stencilWeights
}// namespace hhg