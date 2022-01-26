/*
 * Copyright (c) 2017-2019 Benjamin Mann.
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
#include "hyteg/p1functionspace/variablestencil/VertexDoFVariableStencil.hpp"
#include <hyteg/p2functionspace/variablestencil/P2VariableStencilCommon.hpp>


namespace hyteg {
namespace P2toP1 {

using P2::NumStencilentries2D;

namespace variablestencil {


inline void vertexDoFStencilFromElMat(const Matrixr<3, 6>& elMat, const P2Elements::P2Element& el,
                                      std::array<real_t, NumStencilentries2D::VtV>& VtVStencil,
                                      std::array<real_t, NumStencilentries2D::EtV>& EtVStencil)
{
   for (uint_t i = 0; i < 3; ++i)
   {
      VtVStencil[vertexdof::stencilIndexFromVertex(el[i])] += elMat(0, i);
      EtVStencil[edgedof::stencilIndexFromVertex(el[i + 3])] += elMat(0, i + 3);
   }
}

namespace macroface {

template <class P2ToP1Form>
inline void assembleStencil(const P2ToP1Form& form, const Point3D& x,
                            const Point3D& dirS, const Point3D& dirSE, const Point3D& dirE,
                            const Point3D& dirN, const Point3D& dirNW, const Point3D& dirW,
                            std::array<real_t, NumStencilentries2D::VtV>& VtVStencil,
                            std::array<real_t, NumStencilentries2D::EtV>& EtVStencil)
{
   // nbr face coordinates
   const std::array<Point3D, 3> NE = {x, x + dirE, x + dirN};
   const std::array<Point3D, 3> N = {x, x + dirN, x + dirNW};
   const std::array<Point3D, 3> NW = {x, x + dirNW, x + dirW};
   const std::array<Point3D, 3> SW = {x, x + dirW, x + dirS};
   const std::array<Point3D, 3> S = {x, x + dirS, x + dirSE};
   const std::array<Point3D, 3> SE = {x, x + dirSE, x + dirE};

   VtVStencil.fill(0.0);
   EtVStencil.fill(0.0);

   Matrixr<3,6> elMat;

   // assemble stencils
   form.integrateAll(NE, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementNE_reord, VtVStencil, EtVStencil);

   form.integrateAll(N, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementN_reord, VtVStencil, EtVStencil);

   form.integrateAll(NW, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementNW_reord, VtVStencil, EtVStencil);

   form.integrateAll(SW, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementSW_reord, VtVStencil, EtVStencil);

   form.integrateAll(S, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementS_reord, VtVStencil, EtVStencil);

   form.integrateAll(SE, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementSE_reord, VtVStencil, EtVStencil);

}

template <class P2ToP1Form>
inline void applyVariableStencil(uint_t level,
                                 const Face& face,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Face>& srcVertexDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Face>& srcEdgeDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Face>& dstVertexDoFID,
                                 UpdateType update)
{
   const real_t* srcVertexDoF = face.getData(srcVertexDoFID)->getPointer(level);
   const real_t* srcEdgeDoF   = face.getData(srcEdgeDoFID)->getPointer(level);
   real_t* dstVertexDoF       = face.getData(dstVertexDoFID)->getPointer(level);

   Point3D x0(face.coords[0]), x;
   real_t  h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));

   Point3D d0 = h * (face.coords[1] - face.coords[0]);
   Point3D d2 = h * (face.coords[2] - face.coords[0]);

   P2ToP1Form form;
   form.setGeometryMap(face.getGeometryMap());

   // directions
   const Point3D dirS  = -d2;
   const Point3D dirSE = d0 - d2;
   const Point3D dirE  = d0;
   const Point3D dirW  = -dirE;
   const Point3D dirNW = -dirSE;
   const Point3D dirN  = -dirS;

   // stencil entries
   std::array<real_t, NumStencilentries2D::VtV> vertexToVertexStencil;
   std::array<real_t, NumStencilentries2D::EtV> edgeToVertexStencil;

   // loop over all DOFs
   for (const auto& it : hyteg::edgedof::macroface::Iterator(level, 0))
   {
      x = x0 + walberla::real_c(it.row()) * d2 + walberla::real_c(it.col()) * d0;

      assembleStencil(form, x, dirS, dirSE, dirE, dirN, dirNW, dirW,
                      vertexToVertexStencil, edgeToVertexStencil);

      P2::macroface::applyStencil_DoF(level, it,
                                      vertexToVertexStencil.data(), edgeToVertexStencil.data(),
                                      nullptr, nullptr,
                                      srcVertexDoF, srcEdgeDoF, dstVertexDoF, nullptr, update);
   }
}

} // namespace macroface

namespace macroedge {

template <class P2ToP1Form>
inline void assembleStencil(const P2ToP1Form& formS, const P2ToP1Form& formN, const bool hasNorth, const Point3D& x,
                            const Point3D& dirS, const Point3D& dirSE, const Point3D& dirE,
                            const Point3D& dirN, const Point3D& dirNW, const Point3D& dirW,
                            std::array<real_t, NumStencilentries2D::VtV>& VtVStencil,
                            std::array<real_t, NumStencilentries2D::EtV>& EtVStencil)
{
   // nbr face coordinates
   const std::array<Point3D, 3> NE = {x, x + dirE, x + dirN};
   const std::array<Point3D, 3> N = {x, x + dirN, x + dirNW};
   const std::array<Point3D, 3> NW = {x, x + dirNW, x + dirW};
   const std::array<Point3D, 3> SW = {x, x + dirW, x + dirS};
   const std::array<Point3D, 3> S = {x, x + dirS, x + dirSE};
   const std::array<Point3D, 3> SE = {x, x + dirSE, x + dirE};

   VtVStencil.fill(0.0);
   EtVStencil.fill(0.0);

   Matrixr<3,6> elMat;

   // assemble stencils
   formS.integrateAll(SW, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementSW_reord, VtVStencil, EtVStencil);

   formS.integrateAll(S, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementS_reord, VtVStencil, EtVStencil);

   formS.integrateAll(SE, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementSE_reord, VtVStencil, EtVStencil);

   if (hasNorth)
   {
      formN.integrateAll(NE, elMat);
      vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementNE_reord, VtVStencil, EtVStencil);

      formN.integrateAll(N, elMat);
      vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementN_reord, VtVStencil, EtVStencil);

      formN.integrateAll(NW, elMat);
      vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementNW_reord, VtVStencil, EtVStencil);
   }
}

template <class P2ToP1Form>
inline void applyVariableStencil(uint_t level,
                                 const Edge& edge,
                                 const std::shared_ptr<PrimitiveStorage>& storage,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Edge>& srcVertexDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Edge>& srcEdgeDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Edge>& dstVertexDoFID,
                                 UpdateType update)
{
   typedef stencilDirection SD;

   const real_t* srcVertexDoF = edge.getData(srcVertexDoFID)->getPointer(level);
   const real_t* srcEdgeDoF   = edge.getData(srcEdgeDoFID)->getPointer(level);
   real_t* dstVertexDoF       = edge.getData(dstVertexDoFID)->getPointer(level);

   real_t  h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));

   // faces
   Face* faceS = storage->getFace(edge.neighborFaces()[0]);
   Face* faceN = (edge.getNumNeighborFaces() == 2) ? storage->getFace(edge.neighborFaces()[1]) : nullptr;

   uint_t s_south = 0, e_south = 0, o_south = 0, s_north = 0, e_north = 0, o_north = 0;
   s_south = faceS->vertex_index(edge.neighborVertices()[0]);
   e_south = faceS->vertex_index(edge.neighborVertices()[1]);
   o_south = faceS->vertex_index(faceS->get_vertex_opposite_to_edge(edge.getID()));

   if (faceN)
   {
      s_north = faceN->vertex_index(edge.neighborVertices()[0]);
      e_north = faceN->vertex_index(edge.neighborVertices()[1]);
      o_north = faceN->vertex_index(faceN->get_vertex_opposite_to_edge(edge.getID()));
   }

   Point3D dS_se = h * (faceS->coords[e_south] - faceS->coords[s_south]);
   Point3D dS_so = h * (faceS->coords[o_south] - faceS->coords[s_south]);
   Point3D dS_oe = h * (faceS->coords[e_south] - faceS->coords[o_south]);

   // directions
   Point3D dirS  = -dS_oe;
   Point3D dirE  = dS_se;
   Point3D dirSE = dS_so;
   Point3D dirW  = -dS_se;
   Point3D dirN;
   Point3D dirNW;

   if (faceN)
   {
      Point3D dN_so = h * (faceN->coords[o_north] - faceN->coords[s_north]);
      Point3D dN_oe = h * (faceN->coords[e_north] - faceN->coords[o_north]);
      dirN  = dN_so;
      dirNW = -dN_oe;
   }

   // coords
   Point3D x0  = edge.getCoordinates()[0], x;
   Point3D dx = h * edge.getDirection();

   // forms
   P2ToP1Form formN, formS;
   formS.setGeometryMap(faceS->getGeometryMap());

   if (faceN) formN.setGeometryMap(faceN->getGeometryMap());

   real_t tmp = walberla::real_c(0);

   // stencil entries
   std::array<real_t, NumStencilentries2D::VtV> vertexToVertexStencil;
   std::array<real_t, NumStencilentries2D::EtV> edgeToVertexStencil;

   // loop over all DOFs
   for (const auto& it : hyteg::edgedof::macroedge::Iterator(level, 0))
   {
      uint_t i = it.col();
      x = x0 + walberla::real_c(i) * dx;
      assembleStencil(formS, formN, bool(faceN), x, dirS, dirSE, dirE, dirN, dirNW, dirW, vertexToVertexStencil, edgeToVertexStencil);

      // VERTEX DoF
      if (i != 0)
      {
         if (update == Replace)
         {
            tmp = walberla::real_c(0);
         }
         else
         {
            tmp = dstVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, SD::VERTEX_C)];
         }

         /// on edge vertex dof
         tmp += srcVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, SD::VERTEX_C)] * vertexToVertexStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_C)];

         for (const auto& dir : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF)
         {
            tmp += srcVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, dir)] * vertexToVertexStencil[vertexdof::stencilIndexFromVertex(dir)];
         }

         /// on edge edge dof
         for (const auto& dir : edgedof::macroedge::neighborsOnEdgeFromVertex)
         {
            tmp += srcEdgeDoF[edgedof::macroedge::indexFromVertex(level, i, dir)] * edgeToVertexStencil[edgedof::stencilIndexFromVertex(dir)];
         }

         /// south face vertex dof
         for (const auto& dir : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF)
         {
            tmp += srcVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, dir)] * vertexToVertexStencil[vertexdof::stencilIndexFromVertex(dir)];
         }

         /// south face edge dof
         for (const auto& dir : edgedof::macroedge::neighborsOnSouthFaceFromVertex)
         {
            tmp += srcEdgeDoF[edgedof::macroedge::indexFromVertex(level, i, dir)] * edgeToVertexStencil[edgedof::stencilIndexFromVertex(dir)];
         }

         if (faceN)
         {
            /// north face vertex dof
            for (const auto& dir : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF)
            {
               tmp += srcVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, dir)] * vertexToVertexStencil[vertexdof::stencilIndexFromVertex(dir)];
            }

            /// north face edge dof
            for (const auto& dir : edgedof::macroedge::neighborsOnNorthFaceFromVertex)
            {
               tmp += srcEdgeDoF[edgedof::macroedge::indexFromVertex(level, i, dir)] * edgeToVertexStencil[edgedof::stencilIndexFromVertex(dir)];
            }
         }

         dstVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, SD::VERTEX_C)] = tmp;
      }
   }
}

} // namespace macroedge

namespace macrovertex {

template <class P2ToP1Form>
inline void assembleStencil(const uint_t level, const Vertex& vertex, const std::shared_ptr<PrimitiveStorage>& storage, std::vector<real_t>& VtVStencil, std::vector<real_t>& EtVStencil)
{
   P2ToP1Form form;
   Point3D x, d0, d2;
   Matrixr<3,6> elMat;
   const real_t  h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));
   const uint_t NE = vertex.getNumNeighborEdges();

   std::fill(VtVStencil.begin(), VtVStencil.end(), 0.0);
   std::fill(EtVStencil.begin(), EtVStencil.end(), 0.0);

   for (auto& faceId : vertex.neighborFaces())
   {
      Face* face = storage->getFace(faceId);
      uint_t face_idx = NE + vertex.face_index(face->getID());
      form.setGeometryMap(face->getGeometryMap());

      uint_t vtx = face->vertex_index(vertex.getID());
      auto adj_edges = face->adjacent_edges(vertex.getID());
      uint_t ne = adj_edges.size();

      x = face->coords[vtx];
      d0 = (face->coords[face->vertex_index(storage->getEdge(adj_edges[0])->get_opposite_vertex(vertex.getID()))] - x) * h;
      d2 = (face->coords[face->vertex_index(storage->getEdge(adj_edges[1])->get_opposite_vertex(vertex.getID()))] - x) * h;

      form.integrateAll({{x, x + d0, x + d2}}, elMat);

      // center contribution
      VtVStencil[0] += elMat(0, 0);

      // neighbour edge contribution
      for (uint_t i = 0; i < ne; ++i)
      {
         uint_t edge_idx = vertex.edge_index(adj_edges[i]);
         EtVStencil[edge_idx] += elMat(0, 5 - i);
         VtVStencil[edge_idx + 1] += elMat(0, i + 1);
      }

      // opposite edge contribution
      EtVStencil[face_idx] += elMat(0, 3);
   }

}

template <class P2ToP1Form>
inline void applyVariableStencil(uint_t level,
                                 const Vertex& vertex,
                                 const std::shared_ptr<PrimitiveStorage>& storage,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Vertex>& srcVertexDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Vertex>& srcEdgeDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Vertex>& dstVertexDoFID,
                                 UpdateType update)
{
   const real_t* srcVertexDoF = vertex.getData(srcVertexDoFID)->getPointer(level);
   const real_t* srcEdgeDoF   = vertex.getData(srcEdgeDoFID)->getPointer(level);
   real_t* dstVertexDoF       = vertex.getData(dstVertexDoFID)->getPointer(level);

   const uint_t NE = vertex.getNumNeighborEdges();
   const uint_t NF = vertex.getNumNeighborFaces();

   std::vector<real_t> vertexToVertexStencil(NE + 1);
   std::vector<real_t> edgeToVertexStencil(NE + NF);
   assembleStencil<P2ToP1Form>(level, vertex, storage, vertexToVertexStencil, edgeToVertexStencil);

   if (update == Replace)
   {
      dstVertexDoF[0] = walberla::real_c(0);
   }

   for (size_t i = 0; i < NE + 1; ++i)
   {
      dstVertexDoF[0] += vertexToVertexStencil[i] * srcVertexDoF[i];
   }

   for (size_t i = 0; i < NE + NF; ++i)
   {
      dstVertexDoF[0] += edgeToVertexStencil[i] * srcEdgeDoF[i];
   }
}

} // namespace macrovertex
} // namespace variablestencil
} // namespace P2toP1


namespace P1toP2 {

using P2::NumStencilentries2D;

namespace variablestencil {

// todo
namespace edgeDoFstencilIndices {
/// vertex-edge
const std::array<uint_t, 3> HO_v_N = {vertexdof::stencilIndexFromHorizontalEdge(stencilDirection::VERTEX_W),
                                      vertexdof::stencilIndexFromHorizontalEdge(stencilDirection::VERTEX_E),
                                      vertexdof::stencilIndexFromHorizontalEdge(stencilDirection::VERTEX_NW)
                                     };
const std::array<uint_t, 3> HO_v_S = {vertexdof::stencilIndexFromHorizontalEdge(stencilDirection::VERTEX_W),
                                      vertexdof::stencilIndexFromHorizontalEdge(stencilDirection::VERTEX_E),
                                      vertexdof::stencilIndexFromHorizontalEdge(stencilDirection::VERTEX_SE)
                                     };
const std::array<uint_t, 3> VE_v_W = {vertexdof::stencilIndexFromVerticalEdge(stencilDirection::VERTEX_N),
                                      vertexdof::stencilIndexFromVerticalEdge(stencilDirection::VERTEX_S),
                                      vertexdof::stencilIndexFromVerticalEdge(stencilDirection::VERTEX_NW)
                                     };
const std::array<uint_t, 3> VE_v_E = {vertexdof::stencilIndexFromVerticalEdge(stencilDirection::VERTEX_N),
                                      vertexdof::stencilIndexFromVerticalEdge(stencilDirection::VERTEX_S),
                                      vertexdof::stencilIndexFromVerticalEdge(stencilDirection::VERTEX_SE)
                                     };
const std::array<uint_t, 3> DI_v_SW = {vertexdof::stencilIndexFromDiagonalEdge(stencilDirection::VERTEX_NW),
                                       vertexdof::stencilIndexFromDiagonalEdge(stencilDirection::VERTEX_SE),
                                       vertexdof::stencilIndexFromDiagonalEdge(stencilDirection::VERTEX_SW)
                                      };
const std::array<uint_t, 3> DI_v_NE = {vertexdof::stencilIndexFromDiagonalEdge(stencilDirection::VERTEX_NW),
                                       vertexdof::stencilIndexFromDiagonalEdge(stencilDirection::VERTEX_SE),
                                       vertexdof::stencilIndexFromDiagonalEdge(stencilDirection::VERTEX_NE)
                                      };
// /// edge-edge
// const std::array<uint_t, 3> HO_e_N = {edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_DI_N),
//                                       edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_VE_NW),
//                                       edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_HO_C)
//                                      };
// const std::array<uint_t, 3> HO_e_S = {edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_VE_SE),
//                                       edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_DI_S),
//                                       edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_HO_C)
//                                      };
// const std::array<uint_t, 3> VE_e_W = {edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_DI_W),
//                                       edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_HO_NW),
//                                       edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_VE_C)
//                                      };
// const std::array<uint_t, 3> VE_e_E = {edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_HO_SE),
//                                       edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_DI_E),
//                                       edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_VE_C)
//                                      };
// const std::array<uint_t, 3> DI_e_SW = {edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_HO_S),
//                                        edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_VE_W),
//                                        edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_DI_C)
//                                       };
// const std::array<uint_t, 3> DI_e_NE = {edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_VE_E),
//                                        edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_HO_N),
//                                        edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_DI_C)
//                                       };
}

inline void vertexDoFStencilFromElMat(const Matrixr<6,3>& elMat, const P2Elements::P2Element& el,
                                      std::array<real_t, NumStencilentries2D::VtV>& VtVStencil)
{
   for (uint_t i = 0; i < 3; ++i)
   {
      VtVStencil[vertexdof::stencilIndexFromVertex(el[i])] += elMat(0, i);
   }
}

inline void edgeDoFStencilFromElMat(const Matrixr<6,3>& elMat,
                                    const std::array<uint_t, 3>& idxVtE,
                                    std::array<real_t, NumStencilentries2D::VtE>& VtEStencil)
{
   for (uint_t i = 0; i < 3; ++i)
   {
      VtEStencil[idxVtE[i]] += elMat(5, i);
   }
}

namespace macroface {

template <class P1ToP2Form>
inline void assembleStencil(const P1ToP2Form& form, const Point3D& x,
                            const Point3D& dirS, const Point3D& dirSE, const Point3D& dirE,
                            const Point3D& dirN, const Point3D& dirNW, const Point3D& dirW, const Point3D& dirNE,
                            std::array<real_t, NumStencilentries2D::VtV>& VtVStencil,
                            std::array<real_t, NumStencilentries2D::VtE>& VtEStencil)
{
   // nbr face coordinates
   const std::array<Point3D, 3> NE = {x, x + dirE, x + dirN};
   const std::array<Point3D, 3> N = {x, x + dirN, x + dirNW};
   const std::array<Point3D, 3> NW = {x, x + dirNW, x + dirW};
   const std::array<Point3D, 3> SW = {x, x + dirW, x + dirS};
   const std::array<Point3D, 3> S = {x, x + dirS, x + dirSE};
   const std::array<Point3D, 3> SE = {x, x + dirSE, x + dirE};
   // const std::array<Point3D,3> HO_N = NE;
   const std::array<Point3D, 3> HO_S = {x, x + dirE, x + dirSE};
   const std::array<Point3D, 3> VE_E = {x + dirN, x, x + dirE};
   const std::array<Point3D, 3> VE_W = {x + dirN, x, x + dirNW};
   const std::array<Point3D, 3> DI_SW = {x + dirN, x + dirE, x};
   const std::array<Point3D, 3> DI_NE = {x + dirN, x + dirE, x + dirNE};

   VtVStencil.fill(0.0);
   VtEStencil.fill(0.0);

   Matrixr<6,3> elMat;

   // assemble stencils
   form.integrateAll(NE, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementNE_reord, VtVStencil);
   edgeDoFStencilFromElMat(elMat, edgeDoFstencilIndices::HO_v_N, VtEStencil);

   form.integrateAll(N, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementN_reord, VtVStencil);

   form.integrateAll(NW, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementNW_reord, VtVStencil);

   form.integrateAll(SW, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementSW_reord, VtVStencil);

   form.integrateAll(S, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementS_reord, VtVStencil);

   form.integrateAll(SE, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementSE_reord, VtVStencil);

   form.integrateAll(HO_S, elMat);
   edgeDoFStencilFromElMat(elMat, edgeDoFstencilIndices::HO_v_S, VtEStencil);

   form.integrateAll(VE_E, elMat);
   edgeDoFStencilFromElMat(elMat, edgeDoFstencilIndices::VE_v_E, VtEStencil);

   form.integrateAll(VE_W, elMat);
   edgeDoFStencilFromElMat(elMat, edgeDoFstencilIndices::VE_v_W, VtEStencil);

   form.integrateAll(DI_SW, elMat);
   edgeDoFStencilFromElMat(elMat, edgeDoFstencilIndices::DI_v_SW, VtEStencil);

   form.integrateAll(DI_NE, elMat);
   edgeDoFStencilFromElMat(elMat, edgeDoFstencilIndices::DI_v_NE, VtEStencil);
}

template <class P1ToP2Form>
inline void applyVariableStencil(uint_t level,
                                 const Face& face,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Face>& srcVertexDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Face>& dstVertexDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Face>& dstEdgeDoFID,
                                 UpdateType update)
{
   const real_t* srcVertexDoF = face.getData(srcVertexDoFID)->getPointer(level);
   real_t* dstVertexDoF       = face.getData(dstVertexDoFID)->getPointer(level);
   real_t* dstEdgeDoF         = face.getData(dstEdgeDoFID)->getPointer(level);

   Point3D x0(face.coords[0]), x;
   real_t  h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));

   Point3D d0 = h * (face.coords[1] - face.coords[0]);
   Point3D d2 = h * (face.coords[2] - face.coords[0]);

   P1ToP2Form form;
   form.setGeometryMap(face.getGeometryMap());

   // directions
   const Point3D dirS  = -d2;
   const Point3D dirSE = d0 - d2;
   const Point3D dirE  = d0;
   const Point3D dirW  = -dirE;
   const Point3D dirNW = -dirSE;
   const Point3D dirN  = -dirS;
   const Point3D dirNE = dirN + dirE;

   // stencil entries
   std::array<real_t, NumStencilentries2D::VtV> vertexToVertexStencil;
   std::array<real_t, NumStencilentries2D::VtE> vertexToEdgeStencil;

   // loop over all DOFs
   for (const auto& it : hyteg::edgedof::macroface::Iterator(level, 0))
   {
      x = x0 + walberla::real_c(it.row()) * d2 + walberla::real_c(it.col()) * d0;

      assembleStencil(form, x, dirS, dirSE, dirE, dirN, dirNW, dirW, dirNE,
                      vertexToVertexStencil, vertexToEdgeStencil);

      P2::macroface::applyStencil_DoF(level, it,
                                      vertexToVertexStencil.data(), nullptr,
                                      vertexToEdgeStencil.data(), nullptr,
                                      srcVertexDoF, nullptr, dstVertexDoF, dstEdgeDoF, update);
   }
}

} // namespace macroface

namespace macroedge {

template <class P1ToP2Form>
inline void assembleStencil(const P1ToP2Form& formS, const P1ToP2Form& formN, const bool hasNorth, const Point3D& x,
                            const Point3D& dirS, const Point3D& dirSE, const Point3D& dirE,
                            const Point3D& dirN, const Point3D& dirNW, const Point3D& dirW,
                            std::array<real_t, NumStencilentries2D::VtV>& VtVStencil,
                            std::array<real_t, NumStencilentries2D::VtE>& VtEStencil)
{
   // nbr face coordinates
   const std::array<Point3D, 3> NE = {x, x + dirE, x + dirN};
   const std::array<Point3D, 3> N = {x, x + dirN, x + dirNW};
   const std::array<Point3D, 3> NW = {x, x + dirNW, x + dirW};
   const std::array<Point3D, 3> SW = {x, x + dirW, x + dirS};
   const std::array<Point3D, 3> S = {x, x + dirS, x + dirSE};
   const std::array<Point3D, 3> SE = {x, x + dirSE, x + dirE};
   // const std::array<Point3D,3> HO_N = NE;
   const std::array<Point3D, 3> HO_S = {x, x + dirE, x + dirSE};

   VtVStencil.fill(0.0);
   VtEStencil.fill(0.0);

   Matrixr<6,3> elMat;

   // assemble stencils
   formS.integrateAll(SW, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementSW_reord, VtVStencil);

   formS.integrateAll(S, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementS_reord, VtVStencil);

   formS.integrateAll(SE, elMat);
   vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementSE_reord, VtVStencil);

   formS.integrateAll(HO_S, elMat);
   edgeDoFStencilFromElMat(elMat, edgeDoFstencilIndices::HO_v_S, VtEStencil);

   if (hasNorth)
   {
      formN.integrateAll(NE, elMat);
      vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementNE_reord, VtVStencil);
      edgeDoFStencilFromElMat(elMat, edgeDoFstencilIndices::HO_v_N, VtEStencil);

      formN.integrateAll(N, elMat);
      vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementN_reord, VtVStencil);

      formN.integrateAll(NW, elMat);
      vertexDoFStencilFromElMat(elMat, P2Elements::P2Face::elementNW_reord, VtVStencil);
   }
}

template <class P1ToP2Form>
inline void applyVariableStencil(uint_t level,
                                 const Edge& edge,
                                 const std::shared_ptr<PrimitiveStorage>& storage,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Edge>& srcVertexDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Edge>& dstVertexDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Edge>& dstEdgeDoFID,
                                 UpdateType update)
{
   typedef stencilDirection SD;

   const real_t* srcVertexDoF = edge.getData(srcVertexDoFID)->getPointer(level);
   real_t* dstVertexDoF       = edge.getData(dstVertexDoFID)->getPointer(level);
   real_t* dstEdgeDoF         = edge.getData(dstEdgeDoFID)->getPointer(level);

   real_t  h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));

   // faces
   Face* faceS = storage->getFace(edge.neighborFaces()[0]);
   Face* faceN = (edge.getNumNeighborFaces() == 2) ? storage->getFace(edge.neighborFaces()[1]) : nullptr;

   uint_t s_south = 0, e_south = 0, o_south = 0, s_north = 0, e_north = 0, o_north = 0;
   s_south = faceS->vertex_index(edge.neighborVertices()[0]);
   e_south = faceS->vertex_index(edge.neighborVertices()[1]);
   o_south = faceS->vertex_index(faceS->get_vertex_opposite_to_edge(edge.getID()));

   if (faceN)
   {
      s_north = faceN->vertex_index(edge.neighborVertices()[0]);
      e_north = faceN->vertex_index(edge.neighborVertices()[1]);
      o_north = faceN->vertex_index(faceN->get_vertex_opposite_to_edge(edge.getID()));
   }

   Point3D dS_se = h * (faceS->coords[e_south] - faceS->coords[s_south]);
   Point3D dS_so = h * (faceS->coords[o_south] - faceS->coords[s_south]);
   Point3D dS_oe = h * (faceS->coords[e_south] - faceS->coords[o_south]);

   // directions
   Point3D dirS  = -dS_oe;
   Point3D dirE  = dS_se;
   Point3D dirSE = dS_so;
   Point3D dirW  = -dS_se;
   Point3D dirN;
   Point3D dirNW;

   if (faceN)
   {
      Point3D dN_so = h * (faceN->coords[o_north] - faceN->coords[s_north]);
      Point3D dN_oe = h * (faceN->coords[e_north] - faceN->coords[o_north]);
      dirN  = dN_so;
      dirNW = -dN_oe;
   }

   // coords
   Point3D x0  = edge.getCoordinates()[0], x;
   Point3D dx = h * edge.getDirection();

   // forms
   P1ToP2Form formN, formS;
   formS.setGeometryMap(faceS->getGeometryMap());

   if (faceN) formN.setGeometryMap(faceN->getGeometryMap());

   real_t tmp = walberla::real_c(0);

   // stencil entries
   std::array<real_t, NumStencilentries2D::VtV> vertexToVertexStencil;
   std::array<real_t, NumStencilentries2D::VtE> vertexToEdgeStencil;

   // loop over all DOFs
   for (const auto& it : hyteg::edgedof::macroedge::Iterator(level, 0))
   {
      uint_t i = it.col();
      x = x0 + walberla::real_c(i) * dx;
      assembleStencil(formS, formN, bool(faceN), x, dirS, dirSE, dirE, dirN, dirNW, dirW, vertexToVertexStencil, vertexToEdgeStencil);

      // VERTEX DoF
      if (i != 0)
      {
         if (update == Replace)
         {
            tmp = walberla::real_c(0);
         }
         else
         {
            tmp = dstVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, SD::VERTEX_C)];
         }

         /// on edge vertex dof
         tmp += srcVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, SD::VERTEX_C)] * vertexToVertexStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_C)];

         for (const auto& dir : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF)
         {
            tmp += srcVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, dir)] * vertexToVertexStencil[vertexdof::stencilIndexFromVertex(dir)];
         }

         /// south face vertex dof
         for (const auto& dir : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF)
         {
            tmp += srcVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, dir)] * vertexToVertexStencil[vertexdof::stencilIndexFromVertex(dir)];
         }

         if (faceN)
         {
            /// north face vertex dof
            for (const auto& dir : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF)
            {
               tmp += srcVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, dir)] * vertexToVertexStencil[vertexdof::stencilIndexFromVertex(dir)];
            }
         }

         dstVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, SD::VERTEX_C)] = tmp;
      }

      ////////// (HORIZONTAL) EDGE //////////

      if (update == Replace)
      {
         tmp = walberla::real_c(0);
      }
      else
      {
         tmp = dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge(level, i, SD::EDGE_HO_C)];
      }

      /// on edge
      for (const auto& dir : vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF)
      {
         tmp += srcVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge(level, i, dir)] * vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge(dir)];
      }

      /// on south face
      for (const auto& dir : vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF)
      {
         tmp += srcVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge(level, i, dir)] * vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge(dir)];
      }

      /// on north face
      if (faceN)
      {
         for (const auto& dir : vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF)
         {
            tmp += srcVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge(level, i, dir)] * vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge(dir)];
         }
      }

      dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge(level, i, SD::EDGE_HO_C)] = tmp;
   }
}

} // namespace macroedge

namespace macrovertex {

template <class P1ToP2Form>
inline void assembleStencil(const uint_t level, const Vertex& vertex, const std::shared_ptr<PrimitiveStorage>& storage, std::vector<real_t>& VtVStencil)
{
   P1ToP2Form form;
   Point3D x, d0, d2;
   Matrixr<6,3> elMat;
   const real_t  h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));
   const uint_t NE = vertex.getNumNeighborEdges();

   std::fill(VtVStencil.begin(), VtVStencil.end(), 0.0);

   for (auto& faceId : vertex.neighborFaces())
   {
      Face* face = storage->getFace(faceId);
      uint_t face_idx = NE + vertex.face_index(face->getID());
      form.setGeometryMap(face->getGeometryMap());

      uint_t vtx = face->vertex_index(vertex.getID());
      auto adj_edges = face->adjacent_edges(vertex.getID());
      uint_t ne = adj_edges.size();

      x = face->coords[vtx];
      d0 = (face->coords[face->vertex_index(storage->getEdge(adj_edges[0])->get_opposite_vertex(vertex.getID()))] - x) * h;
      d2 = (face->coords[face->vertex_index(storage->getEdge(adj_edges[1])->get_opposite_vertex(vertex.getID()))] - x) * h;

      form.integrateAll({{x, x + d0, x + d2}}, elMat);

      // center contribution
      VtVStencil[0] += elMat(0, 0);

      // neighbour edge contribution
      for (uint_t i = 0; i < ne; ++i)
      {
         uint_t edge_idx = vertex.edge_index(adj_edges[i]);
         VtVStencil[edge_idx + 1] += elMat(0, i + 1);
      }
   }

}

template <class P1ToP2Form>
inline void applyVariableStencil(uint_t level,
                                 const Vertex& vertex,
                                 const std::shared_ptr<PrimitiveStorage>& storage,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Vertex>& srcVertexDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Vertex>& dstVertexDoFID,
                                 UpdateType update)
{
   const real_t* srcVertexDoF = vertex.getData(srcVertexDoFID)->getPointer(level);
   real_t* dstVertexDoF       = vertex.getData(dstVertexDoFID)->getPointer(level);

   const uint_t NE = vertex.getNumNeighborEdges();
   const uint_t NF = vertex.getNumNeighborFaces();

   std::vector<real_t> vertexToVertexStencil(NE + 1);
   assembleStencil<P1ToP2Form>(level, vertex, storage, vertexToVertexStencil);

   if (update == Replace)
   {
      dstVertexDoF[0] = walberla::real_c(0);
   }

   for (size_t i = 0; i < NE + 1; ++i)
   {
      dstVertexDoF[0] += vertexToVertexStencil[i] * srcVertexDoF[i];
   }
}

} // namespace macrovertex
} // namespace variablestencil
} // namespace P1toP2

} // namespace hyteg