/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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

namespace hyteg {
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
    vertexToEdgeStencil[10] = - (1.0 + 1.0/3.0);
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
}// namespace hyteg