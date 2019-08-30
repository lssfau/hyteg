#pragma once

#include "tinyhhg_core/p2functionspace/P2Elements.hpp"

namespace hyteg {
namespace P2 {
namespace variablestencil {

//template<class P2Form>
//inline void assembleVertexToVertexStencil(const P2Form &form,
//                                          const std::array<Point3D, 3> &coords,
//                                          const P2Elements::P2Element &directions,
//                                          real_t* opr_data) {
//   Point6D matrixRow;
//
//   form.integrateVertexDoF(coords, matrixRow);
//
//   opr_data[vertexdof::stencilIndexFromVertex(directions[0])] += matrixRow[0];
//   opr_data[vertexdof::stencilIndexFromVertex(directions[1])] += matrixRow[1];
//   opr_data[vertexdof::stencilIndexFromVertex(directions[2])] += matrixRow[2];
//}

template<class P2FormT>
inline void assembleEdgeToVertexStencil(const P2FormT &form,
                                        const std::array<Point3D, 3> &coords,
                                        const P2Elements::P2Element &directions,
                                        real_t* opr_data) {
   Point3D matrixRow;

   form.integrateEdgeToVertex(coords, matrixRow);

   opr_data[edgedof::stencilIndexFromVertex(directions[3])] += matrixRow[0];
   opr_data[edgedof::stencilIndexFromVertex(directions[4])] += matrixRow[1];
   opr_data[edgedof::stencilIndexFromVertex(directions[5])] += matrixRow[2];
}

template<class P2Form>
inline void assembleVertexToEdgeStencil(const P2Form &form,
                                        const std::array<Point3D, 3> &coords,
                                        const std::array<uint_t, 3> stencilIndices,
                                        real_t* opr_data) {
   Point3D matrixRow;

   form.integrateVertexToEdge(coords, matrixRow);

   opr_data[stencilIndices[0]] += matrixRow[0];
   opr_data[stencilIndices[1]] += matrixRow[1];
   opr_data[stencilIndices[2]] += matrixRow[2];
}

template<class P2Form>
inline void assembleEdgeToEdgeStencil(const P2Form &form,
                                      const std::array<Point3D, 3> &coords,
                                      const std::array<uint_t, 3> stencilIndices,
                                      real_t* opr_data) {
   Point3D matrixRow;

   form.integrateEdgeToEdge(coords, matrixRow);

   opr_data[stencilIndices[0]] += matrixRow[0];
   opr_data[stencilIndices[1]] += matrixRow[1];
   opr_data[stencilIndices[2]] += matrixRow[2];
}

}
}
}