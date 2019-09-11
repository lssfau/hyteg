/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes.
 *
 * This file is part of HyTeG
 * (see https://i10git.cs.fau.de/hyteg/hyteg).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include "hyteg/p2functionspace/P2Elements.hpp"

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