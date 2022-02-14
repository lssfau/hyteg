/*
 * Copyright (c) 2020 Benjamin Mann.
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

#include "hyteg/primitives/all.hpp"


namespace hyteg {

using walberla::real_t;

namespace stencil {

// container class for 2D stencil directions
struct Directions2D
{
   Directions2D() {}

   /* initialize directions of neighboring nodes
      @param h             local mesh width
      @param edge          a macro-edge
      @param faceS         pointer to the face south of the edge
      @param faceN         pointer to the face north of the edge
   */
   Directions2D(const real_t h, const Edge& edge, const Face* faceS, const Face* faceN)
   {
      // south face

      uint_t s_south = faceS->vertex_index(edge.neighborVertices()[0]);
      uint_t e_south = faceS->vertex_index(edge.neighborVertices()[1]);
      uint_t o_south = faceS->vertex_index(faceS->get_vertex_opposite_to_edge(edge.getID()));

      Point3D dS_se = h * (faceS->getCoordinates()[e_south] - faceS->getCoordinates()[s_south]);
      Point3D dS_so = h * (faceS->getCoordinates()[o_south] - faceS->getCoordinates()[s_south]);
      Point3D dS_oe = h * (faceS->getCoordinates()[e_south] - faceS->getCoordinates()[o_south]);

      S  = -1.0 * dS_oe;
      E  = dS_se;
      SE = dS_so;
      W  = -1.0 * dS_se;

      // north face

      uint_t  s_north, e_north, o_north;

      if (edge.getNumNeighborFaces() == 2)
      {
         s_north = faceN->vertex_index(edge.neighborVertices()[0]);
         e_north = faceN->vertex_index(edge.neighborVertices()[1]);
         o_north = faceN->vertex_index(faceN->get_vertex_opposite_to_edge(edge.getID()));

         Point3D dN_so = h * (faceN->getCoordinates()[o_north] - faceN->getCoordinates()[s_north]);
         Point3D dN_oe = h * (faceN->getCoordinates()[e_north] - faceN->getCoordinates()[o_north]);

         N  = dN_so;
         NW = -1.0 * dN_oe;
      }

   }

   /* initialize directions of neighboring nodes
      @param h             local mesh width
      @param face          a macro-face
   */
   Directions2D(const real_t h, const Face& face)
   {
      Point3D d0 = h * (face.getCoordinates()[1] - face.getCoordinates()[0]);
      Point3D d2 = h * (face.getCoordinates()[2] - face.getCoordinates()[0]);

      S  = -1.0 * d2;
      SE = d0 - 1.0 * d2;
      E  = d0;
      W  = -1.0 * d0;
      NW = -1.0 * d0 + d2;
      N  = d2;
   }

   Point3D S;
   Point3D SE;
   Point3D E;
   Point3D W;
   Point3D NW;
   Point3D N;
};

} // stencil

} // namespace hyteg
