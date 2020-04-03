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
#include "hyteg/p1functionspace/variablestencil/VertexDoFVariableStencil.hpp"


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
inline void assembleEdgeToVertexStencil(const P2FormT& form,
                                        const std::array<Point3D, 3>& coords,
                                        const P2Elements::P2Element& directions,
                                        real_t* opr_data)
{
   Point3D matrixRow;

   form.integrateEdgeToVertex(coords, matrixRow);

   opr_data[edgedof::stencilIndexFromVertex(directions[3])] += matrixRow[0];
   opr_data[edgedof::stencilIndexFromVertex(directions[4])] += matrixRow[1];
   opr_data[edgedof::stencilIndexFromVertex(directions[5])] += matrixRow[2];
}

template<class P2Form>
inline void assembleVertexToEdgeStencil(const P2Form& form,
                                        const std::array<Point3D, 3>& coords,
                                        const std::array<uint_t, 3> stencilIndices,
                                        real_t* opr_data)
{
   Point3D matrixRow;

   form.integrateVertexToEdge(coords, matrixRow);

   opr_data[stencilIndices[0]] += matrixRow[0];
   opr_data[stencilIndices[1]] += matrixRow[1];
   opr_data[stencilIndices[2]] += matrixRow[2];
}

template<class P2Form>
inline void assembleEdgeToEdgeStencil(const P2Form& form,
                                      const std::array<Point3D, 3>& coords,
                                      const std::array<uint_t, 3> stencilIndices,
                                      real_t* opr_data)
{
   Point3D matrixRow;

   form.integrateEdgeToEdge(coords, matrixRow);

   opr_data[stencilIndices[0]] += matrixRow[0];
   opr_data[stencilIndices[1]] += matrixRow[1];
   opr_data[stencilIndices[2]] += matrixRow[2];
}

// stencil indices
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
/// edge-edge
const std::array<uint_t, 3> HO_e_N = {edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_DI_N),
                                      edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_VE_NW),
                                      edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_HO_C)
                                     };
const std::array<uint_t, 3> HO_e_S = {edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_VE_SE),
                                      edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_DI_S),
                                      edgedof::stencilIndexFromHorizontalEdge(stencilDirection::EDGE_HO_C)
                                     };
const std::array<uint_t, 3> VE_e_W = {edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_DI_W),
                                      edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_HO_NW),
                                      edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_VE_C)
                                     };
const std::array<uint_t, 3> VE_e_E = {edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_HO_SE),
                                      edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_DI_E),
                                      edgedof::stencilIndexFromVerticalEdge(stencilDirection::EDGE_VE_C)
                                     };
const std::array<uint_t, 3> DI_e_SW = {edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_HO_S),
                                       edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_VE_W),
                                       edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_DI_C)
                                      };
const std::array<uint_t, 3> DI_e_NE = {edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_VE_E),
                                       edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_HO_N),
                                       edgedof::stencilIndexFromDiagonalEdge(stencilDirection::EDGE_DI_C)
                                      };

namespace macroface {

template <class P2Form>
inline void applyVariableStencil(uint_t level,
                                 const Face& face,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Face>& srcVertexDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Face>& srcEdgeDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Face>& dstVertexDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Face>& dstEdgeDoFID,
                                 UpdateType update)
{
   typedef stencilDirection SD;

   real_t* srcVertexDoF = face.getData(srcVertexDoFID)->getPointer(level);
   real_t* srcEdgeDoF   = face.getData(srcEdgeDoFID)->getPointer(level);
   real_t* dstVertexDoF = face.getData(dstVertexDoFID)->getPointer(level);
   real_t* dstEdgeDoF   = face.getData(dstEdgeDoFID)->getPointer(level);

   Point3D x0(face.coords[0]), x;
   real_t  h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));

   Point3D d0 = h * (face.coords[1] - face.coords[0]);
   Point3D d2 = h * (face.coords[2] - face.coords[0]);

   P2Form form;
   form.setGeometryMap(face.getGeometryMap());

   real_t tmp = 0;

   // directions
   /// vertex to vertex
   const Point3D dirS  = -d2;
   const Point3D dirSE = d0 - d2;
   const Point3D dirE  = d0;
   const Point3D dirW  = -dirE;
   const Point3D dirNW = -dirSE;
   const Point3D dirN  = -dirS;
   /// edge to vertex
   const Point3D dirHO_W  = -0.5 * d0;
   const Point3D dirHO_NW = -0.5 * d0 + d2;
   const Point3D dirHO_E  = -dirHO_W;
   const Point3D dirHO_SE = -dirHO_NW;
   const Point3D dirVE_N  = 0.5 * d2;
   const Point3D dirVE_NW = -d0 + 0.5 * d2;
   const Point3D dirVE_S  = -dirVE_N;
   const Point3D dirVE_SE = -dirVE_NW;
   const Point3D dirDI_SE = 0.5 * d0 - 0.5 * d2;
   const Point3D dirDI_NE = 0.5 * d0 + 0.5 * d2;
   const Point3D dirDI_NW = -dirDI_SE;
   const Point3D dirDI_SW = -dirDI_NE;

   // stencil entries
   std::vector<real_t> vertexToVertexStencil(7);
   std::vector<real_t> edgeToVertexStencil(12);
   std::vector<real_t> vertexToEdgeStencil(12);
   std::vector<real_t> edgeToEdgeStencil(15);

   // loop over all DOFs
   for (const auto& it : hyteg::edgedof::macroface::Iterator(level, 0))
   {
      uint_t j = it.row();
      uint_t i = it.col();
      std::fill(vertexToVertexStencil.begin(), vertexToVertexStencil.end(), 0.0);
      std::fill(edgeToVertexStencil.begin(), edgeToVertexStencil.end(), 0.0);
      std::fill(edgeToEdgeStencil.begin(), edgeToEdgeStencil.end(), 0.0);
      std::fill(vertexToEdgeStencil.begin(), vertexToEdgeStencil.end(), 0.0);

      ////////// VERTEX //////////
      if (!vertexdof::macroface::isVertexOnBoundary(level, it))
      {
         // assemble stencils
         x = x0 + walberla::real_c(j) * d2 + walberla::real_c(i) * d0;

         /// vertex to vertex
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(form, {x, x + dirW, x + dirS}, P1Elements::P1Elements2D::elementSW, vertexToVertexStencil.data());
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(form, {x, x + dirS, x + dirSE}, P1Elements::P1Elements2D::elementS, vertexToVertexStencil.data());
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(form, {x, x + dirSE, x + dirE}, P1Elements::P1Elements2D::elementSE, vertexToVertexStencil.data());
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(form, {x, x + dirE, x + dirN}, P1Elements::P1Elements2D::elementNE, vertexToVertexStencil.data());
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(form, {x, x + dirN, x + dirNW}, P1Elements::P1Elements2D::elementN, vertexToVertexStencil.data());
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(form, {x, x + dirNW, x + dirW}, P1Elements::P1Elements2D::elementNW, vertexToVertexStencil.data());
         /// edge to vertex
         assembleEdgeToVertexStencil<P2Form>(form, {x, x + dirW, x + dirS}, P2Elements::P2Face::elementSW_reord, edgeToVertexStencil.data());
         assembleEdgeToVertexStencil<P2Form>(form, {x, x + dirS, x + dirSE}, P2Elements::P2Face::elementS_reord, edgeToVertexStencil.data());
         assembleEdgeToVertexStencil<P2Form>(form, {x, x + dirSE, x + dirE}, P2Elements::P2Face::elementSE_reord, edgeToVertexStencil.data());
         assembleEdgeToVertexStencil<P2Form>(form, {x, x + dirE, x + dirN}, P2Elements::P2Face::elementNE_reord, edgeToVertexStencil.data());
         assembleEdgeToVertexStencil<P2Form>(form, {x, x + dirN, x + dirNW}, P2Elements::P2Face::elementN_reord, edgeToVertexStencil.data());
         assembleEdgeToVertexStencil<P2Form>(form, {x, x + dirNW, x + dirW}, P2Elements::P2Face::elementNW_reord, edgeToVertexStencil.data());

         // apply Operator
         if (update == Replace)
         {
            tmp = walberla::real_c(0);
         }
         else
         {
            tmp = dstVertexDoF[vertexdof::macroface::indexFromVertex(level, i, j, SD::VERTEX_C)];
         }

         /// vertex to vertex
         for (const auto& dir : vertexdof::macroface::neighborsWithCenter)
         {
            tmp += srcVertexDoF[vertexdof::macroface::indexFromVertex(level, i, j, dir)] *
                   vertexToVertexStencil[vertexdof::stencilIndexFromVertex(dir)];
         }

         /// edge to vertex
         for (const auto& dir : edgedof::macroface::neighborsFromVertex)
         {
            tmp += srcEdgeDoF[edgedof::macroface::indexFromVertex(level, i, j, dir)] *
                   edgeToVertexStencil[edgedof::stencilIndexFromVertex(dir)];
         }

        dstVertexDoF[vertexdof::macroface::indexFromVertex(level, i, j, SD::VERTEX_C)] = tmp;
      }

      ////////// HORIZONTAL EDGE //////////
      if (!edgedof::macroface::isHorizontalEdgeOnBoundary(level, it))
      {
         // assemble stencils
         x = x0 + walberla::real_c(j) * d2 + walberla::real_c(i + 0.5) * d0;
         /// vertex to edge
         assembleVertexToEdgeStencil<P2Form>(form, {x + dirHO_W, x + dirHO_E, x + dirHO_NW}, HO_v_N, vertexToEdgeStencil.data());
         assembleVertexToEdgeStencil<P2Form>(form, {x + dirHO_W, x + dirHO_E, x + dirHO_SE}, HO_v_S, vertexToEdgeStencil.data());
         /// edge to edge
         assembleEdgeToEdgeStencil<P2Form>(form, {x + dirHO_W, x + dirHO_E, x + dirHO_NW}, HO_e_N, edgeToEdgeStencil.data());
         assembleEdgeToEdgeStencil<P2Form>(form, {x + dirHO_W, x + dirHO_E, x + dirHO_SE}, HO_e_S, edgeToEdgeStencil.data());

         // apply Operator
         if (update == Replace)
         {
            tmp = walberla::real_c(0);
         }
         else
         {
            tmp = dstEdgeDoF[edgedof::macroface::indexFromHorizontalEdge(level, i, j, SD::EDGE_HO_C)];
         }

         /// vertex to edge
         for (const auto& dir : vertexdof::macroface::neighborsFromHorizontalEdge)
         {
            tmp += srcVertexDoF[vertexdof::macroface::indexFromHorizontalEdge(level, i, j, dir)] *
                   vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge(dir)];
         }

         /// edge to edge
         for (const auto& dir : edgedof::macroface::neighborsFromHorizontalEdge)
         {
            tmp += srcEdgeDoF[edgedof::macroface::indexFromHorizontalEdge(level, i, j, dir)] *
                   edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge(dir)];
         }

         dstEdgeDoF[edgedof::macroface::indexFromHorizontalEdge(level, i, j, SD::EDGE_HO_C)] = tmp;
      }

      ////////// VERTICAL EDGE //////////
      if (!edgedof::macroface::isVerticalEdgeOnBoundary(level, it))
      {
         // assemble stencils
         x = x0 + walberla::real_c(j + 0.5) * d2 + walberla::real_c(i) * d0;
         /// vertex to edge
         assembleVertexToEdgeStencil<P2Form>(form, {x + dirVE_N, x + dirVE_S, x + dirVE_NW}, VE_v_W, vertexToEdgeStencil.data());
         assembleVertexToEdgeStencil<P2Form>(form, {x + dirVE_N, x + dirVE_S, x + dirVE_SE}, VE_v_E, vertexToEdgeStencil.data());
         /// edge to edge
         assembleEdgeToEdgeStencil<P2Form>(form, {x + dirVE_N, x + dirVE_S, x + dirVE_NW}, VE_e_W, edgeToEdgeStencil.data());
         assembleEdgeToEdgeStencil<P2Form>(form, {x + dirVE_N, x + dirVE_S, x + dirVE_SE}, VE_e_E, edgeToEdgeStencil.data());

         // apply Operator
         if (update == Replace)
         {
            tmp = walberla::real_c(0);
         }
         else
         {
            tmp = dstEdgeDoF[edgedof::macroface::indexFromVerticalEdge(level, i, j, SD::EDGE_VE_C)];
         }

         /// vertex to edge
         for (const auto& dir : vertexdof::macroface::neighborsFromVerticalEdge)
         {
            tmp += srcVertexDoF[vertexdof::macroface::indexFromVerticalEdge(level, i, j, dir)] *
                   vertexToEdgeStencil[vertexdof::stencilIndexFromVerticalEdge(dir)];
         }

         /// edge to edge
         for (const auto& dir : edgedof::macroface::neighborsFromVerticalEdge)
         {
            tmp += srcEdgeDoF[edgedof::macroface::indexFromVerticalEdge(level, i, j, dir)] *
                   edgeToEdgeStencil[edgedof::stencilIndexFromVerticalEdge(dir)];
         }

         dstEdgeDoF[edgedof::macroface::indexFromVerticalEdge(level, i, j, SD::EDGE_VE_C)] = tmp;
      }

      ////////// DIAGONAL EDGE //////////
      if (!edgedof::macroface::isDiagonalEdgeOnBoundary(level, it))
      {
         // assemble stencils
         x = x0 + walberla::real_c(j + 0.5) * d2 + walberla::real_c(i + 0.5) * d0;
         /// vertex to edge
         assembleVertexToEdgeStencil<P2Form>(form, {x + dirDI_NW, x + dirDI_SE, x + dirDI_SW}, DI_v_SW, vertexToEdgeStencil.data());
         assembleVertexToEdgeStencil<P2Form>(form, {x + dirDI_NW, x + dirDI_SE, x + dirDI_NE}, DI_v_NE, vertexToEdgeStencil.data());
         /// edge to edge
         assembleEdgeToEdgeStencil<P2Form>(form, {x + dirDI_NW, x + dirDI_SE, x + dirDI_SW}, DI_e_SW, edgeToEdgeStencil.data());
         assembleEdgeToEdgeStencil<P2Form>(form, {x + dirDI_NW, x + dirDI_SE, x + dirDI_NE}, DI_e_NE, edgeToEdgeStencil.data());

         // apply Operator
         if (update == Replace)
         {
            tmp = walberla::real_c(0);
         }
         else
         {
            tmp = dstEdgeDoF[edgedof::macroface::indexFromDiagonalEdge(level, i, j, SD::EDGE_DI_C)];
         }

         /// vertex to edge
         for (const auto& dir : vertexdof::macroface::neighborsFromDiagonalEdge)
         {
            tmp += srcVertexDoF[vertexdof::macroface::indexFromDiagonalEdge(level, i, j, dir)] *
                   vertexToEdgeStencil[vertexdof::stencilIndexFromDiagonalEdge(dir)];
         }

         /// edge to edge
         for (const auto& dir : edgedof::macroface::neighborsFromDiagonalEdge)
         {
            tmp += srcEdgeDoF[edgedof::macroface::indexFromDiagonalEdge(level, i, j, dir)] *
                   edgeToEdgeStencil[edgedof::stencilIndexFromDiagonalEdge(dir)];
         }

         dstEdgeDoF[edgedof::macroface::indexFromDiagonalEdge(level, i, j, SD::EDGE_DI_C)] = tmp;
      }
   }
}

template <class P2Form>
inline void smoothGSVariableStencil(uint_t level,
                                    const Face& face,
                                    const PrimitiveDataID<FunctionMemory<real_t>, Face>& dstVertexDoFID,
                                    const PrimitiveDataID<FunctionMemory<real_t>, Face>& dstEdgeDoFID,
                                    const PrimitiveDataID<FunctionMemory<real_t>, Face>& rhsVertexDoFID,
                                    const PrimitiveDataID<FunctionMemory<real_t>, Face>& rhsEdgeDoFID)
{
   typedef stencilDirection SD;

   real_t* dstVertexDoF = face.getData(dstVertexDoFID)->getPointer(level);
   real_t* dstEdgeDoF   = face.getData(dstEdgeDoFID)->getPointer(level);
   real_t* rhsVertexDoF = face.getData(rhsVertexDoFID)->getPointer(level);
   real_t* rhsEdgeDoF   = face.getData(rhsEdgeDoFID)->getPointer(level);

   // uint_t rowsize       = levelinfo::num_microvertices_per_edge(level);
   // uint_t inner_rowsize = rowsize;

   // auto rhs = face.getData(rhsId)->getPointer(level);
   // auto dst = face.getData(dstId)->getPointer(level);

   Point3D x0(face.coords[0]), x;
   real_t  h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));

   Point3D d0 = h * (face.coords[1] - face.coords[0]);
   Point3D d2 = h * (face.coords[2] - face.coords[0]);

   P2Form form;
   form.setGeometryMap(face.getGeometryMap());

   real_t tmp = 0;

   // directions
   /// vertex to vertex
   const Point3D dirS  = -d2;
   const Point3D dirSE = d0 - d2;
   const Point3D dirE  = d0;
   const Point3D dirW  = -dirE;
   const Point3D dirNW = -dirSE;
   const Point3D dirN  = -dirS;
   /// edge to vertex
   const Point3D dirHO_W  = -0.5 * d0;
   const Point3D dirHO_NW = -0.5 * d0 + d2;
   const Point3D dirHO_E  = -dirHO_W;
   const Point3D dirHO_SE = -dirHO_NW;
   const Point3D dirVE_N  = 0.5 * d2;
   const Point3D dirVE_NW = -d0 + 0.5 * d2;
   const Point3D dirVE_S  = -dirVE_N;
   const Point3D dirVE_SE = -dirVE_NW;
   const Point3D dirDI_SE = 0.5 * d0 - 0.5 * d2;
   const Point3D dirDI_NE = 0.5 * d0 + 0.5 * d2;
   const Point3D dirDI_NW = -dirDI_SE;
   const Point3D dirDI_SW = -dirDI_NE;

   // stencil entries
   std::vector<real_t> vertexToVertexStencil(7);
   std::vector<real_t> edgeToVertexStencil(12);
   std::vector<real_t> vertexToEdgeStencil(12);
   std::vector<real_t> edgeToEdgeStencil(15);

   // loop over all DOFs
   // for (uint_t j = 1; j < rowsize - 2; ++j)
   for (const auto& it : hyteg::edgedof::macroface::Iterator(level, 0))
   {
      uint_t j = it.row();
      uint_t i = it.col();
      std::fill(vertexToVertexStencil.begin(), vertexToVertexStencil.end(), 0.0);
      std::fill(edgeToVertexStencil.begin(), edgeToVertexStencil.end(), 0.0);
      std::fill(edgeToEdgeStencil.begin(), edgeToEdgeStencil.end(), 0.0);
      std::fill(vertexToEdgeStencil.begin(), vertexToEdgeStencil.end(), 0.0);

      // x += walberla::real_c(j) * d2 + d0;

      // for (uint_t i = 1; i < inner_rowsize - 2; ++i)
      // {

      ////////// VERTEX //////////
      if (!vertexdof::macroface::isVertexOnBoundary(level, it))
      {
         // assemble stencils
         x = x0 + walberla::real_c(j) * d2 + walberla::real_c(i) * d0;

         /// vertex to vertex
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(form, {x, x + dirW, x + dirS}, P1Elements::P1Elements2D::elementSW, vertexToVertexStencil.data());
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(form, {x, x + dirS, x + dirSE}, P1Elements::P1Elements2D::elementS, vertexToVertexStencil.data());
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(form, {x, x + dirSE, x + dirE}, P1Elements::P1Elements2D::elementSE, vertexToVertexStencil.data());
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(form, {x, x + dirE, x + dirN}, P1Elements::P1Elements2D::elementNE, vertexToVertexStencil.data());
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(form, {x, x + dirN, x + dirNW}, P1Elements::P1Elements2D::elementN, vertexToVertexStencil.data());
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(form, {x, x + dirNW, x + dirW}, P1Elements::P1Elements2D::elementNW, vertexToVertexStencil.data());
         /// edge to vertex
         assembleEdgeToVertexStencil<P2Form>(form, {x, x + dirW, x + dirS}, P2Elements::P2Face::elementSW_reord, edgeToVertexStencil.data());
         assembleEdgeToVertexStencil<P2Form>(form, {x, x + dirS, x + dirSE}, P2Elements::P2Face::elementS_reord, edgeToVertexStencil.data());
         assembleEdgeToVertexStencil<P2Form>(form, {x, x + dirSE, x + dirE}, P2Elements::P2Face::elementSE_reord, edgeToVertexStencil.data());
         assembleEdgeToVertexStencil<P2Form>(form, {x, x + dirE, x + dirN}, P2Elements::P2Face::elementNE_reord, edgeToVertexStencil.data());
         assembleEdgeToVertexStencil<P2Form>(form, {x, x + dirN, x + dirNW}, P2Elements::P2Face::elementN_reord, edgeToVertexStencil.data());
         assembleEdgeToVertexStencil<P2Form>(form, {x, x + dirNW, x + dirW}, P2Elements::P2Face::elementNW_reord, edgeToVertexStencil.data());

         // apply GS
         tmp = rhsVertexDoF[vertexdof::macroface::indexFromVertex(level, i, j, SD::VERTEX_C)];

         // for (uint_t k = 0; k < vertexdof::macroface::neighborsWithoutCenter.size(); ++k)
         // {
         //    tmp -= vertexToVertexStencil[vertexdof::stencilIndexFromVertex(vertexdof::macroface::neighborsWithoutCenter[k])] *
         //           dst[vertexdof::macroface::indexFromVertex(level, i, j, vertexdof::macroface::neighborsWithoutCenter[k])];
         // }
         /// vertex to vertex
         for (const auto& dir : vertexdof::macroface::neighborsWithoutCenter)
         {
            tmp -= dstVertexDoF[vertexdof::macroface::indexFromVertex(level, i, j, dir)] *
                   vertexToVertexStencil[vertexdof::stencilIndexFromVertex(dir)];
         }

         /// edge to vertex
         for (const auto& dir : edgedof::macroface::neighborsFromVertex)
         {
            tmp -= dstEdgeDoF[edgedof::macroface::indexFromVertex(level, i, j, dir)] *
                   edgeToVertexStencil[edgedof::stencilIndexFromVertex(dir)];
         }

         dstVertexDoF[vertexdof::macroface::indexFromVertex(level, i, j, SD::VERTEX_C)] =
            tmp / vertexToVertexStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_C)];
      }

      ////////// HORIZONTAL EDGE //////////
      if (!edgedof::macroface::isHorizontalEdgeOnBoundary(level, it))
      {
         // assemble stencils
         x = x0 + walberla::real_c(j) * d2 + walberla::real_c(i + 0.5) * d0;
         /// vertex to edge
         assembleVertexToEdgeStencil<P2Form>(form, {x + dirHO_W, x + dirHO_E, x + dirHO_NW}, HO_v_N, vertexToEdgeStencil.data());
         assembleVertexToEdgeStencil<P2Form>(form, {x + dirHO_W, x + dirHO_E, x + dirHO_SE}, HO_v_S, vertexToEdgeStencil.data());
         /// edge to edge
         assembleEdgeToEdgeStencil<P2Form>(form, {x + dirHO_W, x + dirHO_E, x + dirHO_NW}, HO_e_N, edgeToEdgeStencil.data());
         assembleEdgeToEdgeStencil<P2Form>(form, {x + dirHO_W, x + dirHO_E, x + dirHO_SE}, HO_e_S, edgeToEdgeStencil.data());

         // apply GS
         tmp = rhsEdgeDoF[edgedof::macroface::indexFromHorizontalEdge(level, i, j, SD::EDGE_HO_C)];

         /// vertex to edge
         for (const auto& dir : vertexdof::macroface::neighborsFromHorizontalEdge)
         {
            tmp -= dstVertexDoF[vertexdof::macroface::indexFromHorizontalEdge(level, i, j, dir)] *
                   vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge(dir)];
         }

         /// edge to edge
         for (const auto& dir : edgedof::macroface::neighborsFromHorizontalEdgeWithoutCenter)
         {
            tmp -= dstEdgeDoF[edgedof::macroface::indexFromHorizontalEdge(level, i, j, dir)] *
                   edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge(dir)];
         }

         dstEdgeDoF[edgedof::macroface::indexFromHorizontalEdge(level, i, j, SD::EDGE_HO_C)] =
            tmp / edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_HO_C)];
      }

      ////////// VERTICAL EDGE //////////
      if (!edgedof::macroface::isVerticalEdgeOnBoundary(level, it))
      {
         // assemble stencils
         x = x0 + walberla::real_c(j + 0.5) * d2 + walberla::real_c(i) * d0;
         /// vertex to edge
         assembleVertexToEdgeStencil<P2Form>(form, {x + dirVE_N, x + dirVE_S, x + dirVE_NW}, VE_v_W, vertexToEdgeStencil.data());
         assembleVertexToEdgeStencil<P2Form>(form, {x + dirVE_N, x + dirVE_S, x + dirVE_SE}, VE_v_E, vertexToEdgeStencil.data());
         /// edge to edge
         assembleEdgeToEdgeStencil<P2Form>(form, {x + dirVE_N, x + dirVE_S, x + dirVE_NW}, VE_e_W, edgeToEdgeStencil.data());
         assembleEdgeToEdgeStencil<P2Form>(form, {x + dirVE_N, x + dirVE_S, x + dirVE_SE}, VE_e_E, edgeToEdgeStencil.data());

         // apply GS
         tmp = rhsEdgeDoF[edgedof::macroface::indexFromVerticalEdge(level, i, j, SD::EDGE_VE_C)];

         /// vertex to edge
         for (const auto& dir : vertexdof::macroface::neighborsFromVerticalEdge)
         {
            tmp -= dstVertexDoF[vertexdof::macroface::indexFromVerticalEdge(level, i, j, dir)] *
                   vertexToEdgeStencil[vertexdof::stencilIndexFromVerticalEdge(dir)];
         }

         /// edge to edge
         for (const auto& dir : edgedof::macroface::neighborsFromVerticalEdgeWithoutCenter)
         {
            tmp -= dstEdgeDoF[edgedof::macroface::indexFromVerticalEdge(level, i, j, dir)] *
                   edgeToEdgeStencil[edgedof::stencilIndexFromVerticalEdge(dir)];
         }

         dstEdgeDoF[edgedof::macroface::indexFromVerticalEdge(level, i, j, SD::EDGE_VE_C)] =
            tmp / edgeToEdgeStencil[edgedof::stencilIndexFromVerticalEdge(SD::EDGE_VE_C)];
      }

      ////////// DIAGONAL EDGE //////////
      if (!edgedof::macroface::isDiagonalEdgeOnBoundary(level, it))
      {
         // assemble stencils
         x = x0 + walberla::real_c(j + 0.5) * d2 + walberla::real_c(i + 0.5) * d0;
         /// vertex to edge
         assembleVertexToEdgeStencil<P2Form>(form, {x + dirDI_NW, x + dirDI_SE, x + dirDI_SW}, DI_v_SW, vertexToEdgeStencil.data());
         assembleVertexToEdgeStencil<P2Form>(form, {x + dirDI_NW, x + dirDI_SE, x + dirDI_NE}, DI_v_NE, vertexToEdgeStencil.data());
         /// edge to edge
         assembleEdgeToEdgeStencil<P2Form>(form, {x + dirDI_NW, x + dirDI_SE, x + dirDI_SW}, DI_e_SW, edgeToEdgeStencil.data());
         assembleEdgeToEdgeStencil<P2Form>(form, {x + dirDI_NW, x + dirDI_SE, x + dirDI_NE}, DI_e_NE, edgeToEdgeStencil.data());

         // apply GS
         tmp = rhsEdgeDoF[edgedof::macroface::indexFromDiagonalEdge(level, i, j, SD::EDGE_DI_C)];

         /// vertex to edge
         for (const auto& dir : vertexdof::macroface::neighborsFromDiagonalEdge)
         {
            tmp -= dstVertexDoF[vertexdof::macroface::indexFromDiagonalEdge(level, i, j, dir)] *
                   vertexToEdgeStencil[vertexdof::stencilIndexFromDiagonalEdge(dir)];
         }

         /// edge to edge
         for (const auto& dir : edgedof::macroface::neighborsFromDiagonalEdgeWithoutCenter)
         {
            tmp -= dstEdgeDoF[edgedof::macroface::indexFromDiagonalEdge(level, i, j, dir)] *
                   edgeToEdgeStencil[edgedof::stencilIndexFromDiagonalEdge(dir)];
         }

         dstEdgeDoF[edgedof::macroface::indexFromDiagonalEdge(level, i, j, SD::EDGE_DI_C)] =
            tmp / edgeToEdgeStencil[edgedof::stencilIndexFromDiagonalEdge(SD::EDGE_DI_C)];
      }

      // x += d0;
      // }

      // --inner_rowsize;
   }
}

} // namespace macroface

namespace macroedge {

template <class P2Form>
inline void applyVariableStencil(uint_t level,
                                 const Edge& edge,
                                 const std::shared_ptr<PrimitiveStorage>& storage,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Edge>& srcVertexDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Edge>& srcEdgeDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Edge>& dstVertexDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Edge>& dstEdgeDoFID,
                                 UpdateType update)
{
   typedef stencilDirection SD;

   real_t* srcVertexDoF = edge.getData(srcVertexDoFID)->getPointer(level);
   real_t* srcEdgeDoF   = edge.getData(srcEdgeDoFID)->getPointer(level);
   real_t* dstVertexDoF = edge.getData(dstVertexDoFID)->getPointer(level);
   real_t* dstEdgeDoF   = edge.getData(dstEdgeDoFID)->getPointer(level);

   real_t  h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));

   // faces
   Face* faceS = storage->getFace(edge.neighborFaces()[0]);
   Face* faceN = (edge.getNumNeighborFaces() == 2) ? storage->getFace(edge.neighborFaces()[1]) : nullptr;

   uint_t s_south, e_south, o_south, s_north, e_north, o_north;
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
   /// vertex to vertex
   Point3D dir_S  = -dS_oe;
   Point3D dir_E  = dS_se;
   Point3D dir_SE = dS_so;
   Point3D dir_W  = -dS_se;
   Point3D dir_N;
   Point3D dir_NW;
   /// edge to vertex
   Point3D dirHO_W  = -0.5 * dS_se;
   Point3D dirHO_E  = -dirHO_W;
   Point3D dirHO_SE = 0.5 * dS_se - dS_oe;
   Point3D dirHO_NW;

   if (faceN)
   {
      Point3D dN_so = h * (faceN->coords[o_north] - faceN->coords[s_north]);
      Point3D dN_oe = h * (faceN->coords[e_north] - faceN->coords[o_north]);
      dir_N  = dN_so;
      dir_NW = -dN_oe;
      dirHO_NW = -0.5 * dS_se + dN_so;
   }

   // const Point3D dirS  = -d2;
   // const Point3D dirSE = d0 - d2;
   // const Point3D dirE  = d0;
   // const Point3D dirW  = -dirW;
   // const Point3D dirNW = -dirSE;
   // const Point3D dirN  = -dirS;
   // Point3D dirVE_N  = 0.5 * d2;
   // Point3D dirVE_NW = -d0 + 0.5 * d2;
   // Point3D dirVE_S  = -dirVE_N;
   // Point3D dirVE_SE = -dirVE_NW;
   // Point3D dirDI_SE = 0.5 * d0 - 0.5 * d2;
   // Point3D dirDI_NE = 0.5 * d0 + 0.5 * d2;
   // Point3D dirDI_NW = -dirDI_SE;
   // Point3D dirDI_SW = -dirDI_NE;

   // coords
   Point3D x0  = edge.getCoordinates()[0], x;
   Point3D dx = h * edge.getDirection();

   // forms
   P2Form formS, formN;
   formS.setGeometryMap(faceS->getGeometryMap());

   if (faceN) formN.setGeometryMap(faceN->getGeometryMap());

   real_t tmp = walberla::real_c(0);

   // stencil entries
   std::vector<real_t> vertexToVertexStencil(7);
   std::vector<real_t> edgeToVertexStencil(12);
   std::vector<real_t> vertexToEdgeStencil(12);
   std::vector<real_t> edgeToEdgeStencil(15);

   // loop over all DOFs
   for (const auto& it : hyteg::edgedof::macroedge::Iterator(level, 0))
   {
      uint_t i = it.col();
      std::fill(vertexToVertexStencil.begin(), vertexToVertexStencil.end(), 0.0);
      std::fill(edgeToVertexStencil.begin(), edgeToVertexStencil.end(), 0.0);
      std::fill(edgeToEdgeStencil.begin(), edgeToEdgeStencil.end(), 0.0);
      std::fill(vertexToEdgeStencil.begin(), vertexToEdgeStencil.end(), 0.0);

      ////////// VERTEX //////////
      if (i != 0)
      {
         // assemble stencils
         x = x0 + walberla::real_c(i) * dx;
         /// south face
         //// vertex to vertex
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(formS, {x, x + dir_W, x + dir_S}, P1Elements::P1Elements2D::elementSW, vertexToVertexStencil.data());
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(formS, {x, x + dir_S, x + dir_SE}, P1Elements::P1Elements2D::elementS, vertexToVertexStencil.data());
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(formS, {x, x + dir_SE, x + dir_E}, P1Elements::P1Elements2D::elementSE, vertexToVertexStencil.data());
         //// edge to vertex
         assembleEdgeToVertexStencil<P2Form>(formS, {x, x + dir_W, x + dir_S}, P2Elements::P2Face::elementSW_reord, edgeToVertexStencil.data());
         assembleEdgeToVertexStencil<P2Form>(formS, {x, x + dir_S, x + dir_SE}, P2Elements::P2Face::elementS_reord, edgeToVertexStencil.data());
         assembleEdgeToVertexStencil<P2Form>(formS, {x, x + dir_SE, x + dir_E}, P2Elements::P2Face::elementSE_reord, edgeToVertexStencil.data());

         /// north face
         if (faceN)
         {
            //// vertex to vertex
            vertexdof::variablestencil::assembleLocalStencil<P2Form>(formN, {x, x + dir_E, x + dir_N}, P1Elements::P1Elements2D::elementNE, vertexToVertexStencil.data());
            vertexdof::variablestencil::assembleLocalStencil<P2Form>(formN, {x, x + dir_N, x + dir_NW}, P1Elements::P1Elements2D::elementN, vertexToVertexStencil.data());
            vertexdof::variablestencil::assembleLocalStencil<P2Form>(formN, {x, x + dir_NW, x + dir_W}, P1Elements::P1Elements2D::elementNW, vertexToVertexStencil.data());
            //// edge to vertex
            assembleEdgeToVertexStencil<P2Form>(formN, {x, x + dir_E, x + dir_N}, P2Elements::P2Face::elementNE_reord, edgeToVertexStencil.data());
            assembleEdgeToVertexStencil<P2Form>(formN, {x, x + dir_N, x + dir_NW}, P2Elements::P2Face::elementN_reord, edgeToVertexStencil.data());
            assembleEdgeToVertexStencil<P2Form>(formN, {x, x + dir_NW, x + dir_W}, P2Elements::P2Face::elementNW_reord, edgeToVertexStencil.data());
         }

         // apply Operator
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

      ////////// (HORIZONTAL) EDGE //////////
      // assemble stencils
      x = x0 + walberla::real_c(i + 0.5) * dx;
      /// south face
      //// vertex to edge
      assembleVertexToEdgeStencil<P2Form>(formS, {x + dirHO_W, x + dirHO_E, x + dirHO_SE}, HO_v_S, vertexToEdgeStencil.data());
      //// edge to edge
      assembleEdgeToEdgeStencil<P2Form>(formS, {x + dirHO_W, x + dirHO_E, x + dirHO_SE}, HO_e_S, edgeToEdgeStencil.data());

      /// north face
      if (faceN)
      {
         //// vertex to edge
         assembleVertexToEdgeStencil<P2Form>(formN, {x + dirHO_W, x + dirHO_E, x + dirHO_NW}, HO_v_N, vertexToEdgeStencil.data());
         //// edge to edge
         assembleEdgeToEdgeStencil<P2Form>(formN, {x + dirHO_W, x + dirHO_E, x + dirHO_NW}, HO_e_N, edgeToEdgeStencil.data());
      }

      // apply Operator
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

      for (const auto& dir : edgedof::macroedge::neighborsOnEdgeFromHorizontalEdge)
      {
         tmp += srcEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge(level, i, dir)] * edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge(dir)];
      }

      /// on south face
      for (const auto& dir : vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF)
      {
         tmp += srcVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge(level, i, dir)] * vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge(dir)];
      }

      for (const auto& dir : edgedof::macroedge::neighborsOnSouthFaceFromHorizontalEdge)
      {
         tmp += srcEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge(level, i, dir)] * edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge(dir)];
      }

      /// on north face
      if (faceN)
      {
         for (const auto& dir : vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF)
         {
            tmp += srcVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge(level, i, dir)] * vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge(dir)];
         }

         for (const auto& dir : edgedof::macroedge::neighborsOnNorthFaceFromHorizontalEdge)
         {
            tmp += srcEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge(level, i, dir)] * edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge(dir)];
         }
      }

      dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge(level, i, SD::EDGE_HO_C)] = tmp;
   }
}

template <class P2Form>
inline void smoothGSVariableStencil(uint_t level,
                                    const Edge& edge,
                                    const std::shared_ptr<PrimitiveStorage>& storage,
                                    const PrimitiveDataID<FunctionMemory<real_t>, Edge>& dstVertexDoFID,
                                    const PrimitiveDataID<FunctionMemory<real_t>, Edge>& dstEdgeDoFID,
                                    const PrimitiveDataID<FunctionMemory<real_t>, Edge>& rhsVertexDoFID,
                                    const PrimitiveDataID<FunctionMemory<real_t>, Edge>& rhsEdgeDoFID)
{
   typedef stencilDirection SD;

   real_t* dstVertexDoF = edge.getData(dstVertexDoFID)->getPointer(level);
   real_t* dstEdgeDoF   = edge.getData(dstEdgeDoFID)->getPointer(level);
   real_t* rhsVertexDoF = edge.getData(rhsVertexDoFID)->getPointer(level);
   real_t* rhsEdgeDoF   = edge.getData(rhsEdgeDoFID)->getPointer(level);

   real_t  h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));

   // faces
   Face* faceS = storage->getFace(edge.neighborFaces()[0]);
   Face* faceN = (edge.getNumNeighborFaces() == 2) ? storage->getFace(edge.neighborFaces()[1]) : nullptr;

   uint_t s_south, e_south, o_south, s_north, e_north, o_north;
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
   /// vertex to vertex
   Point3D dir_S  = -dS_oe;
   Point3D dir_E  = dS_se;
   Point3D dir_SE = dS_so;
   Point3D dir_W  = -dS_se;
   Point3D dir_N;
   Point3D dir_NW;
   /// edge to vertex
   Point3D dirHO_W  = -0.5 * dS_se;
   Point3D dirHO_E  = -dirHO_W;
   Point3D dirHO_SE = 0.5 * dS_se - dS_oe;
   Point3D dirHO_NW;

   if (faceN)
   {
      Point3D dN_so = h * (faceN->coords[o_north] - faceN->coords[s_north]);
      Point3D dN_oe = h * (faceN->coords[e_north] - faceN->coords[o_north]);
      dir_N  = dN_so;
      dir_NW = -dN_oe;
      dirHO_NW = -0.5 * dS_se + dN_so;
   }

   // const Point3D dirS  = -d2;
   // const Point3D dirSE = d0 - d2;
   // const Point3D dirE  = d0;
   // const Point3D dirW  = -dirW;
   // const Point3D dirNW = -dirSE;
   // const Point3D dirN  = -dirS;
   // Point3D dirVE_N  = 0.5 * d2;
   // Point3D dirVE_NW = -d0 + 0.5 * d2;
   // Point3D dirVE_S  = -dirVE_N;
   // Point3D dirVE_SE = -dirVE_NW;
   // Point3D dirDI_SE = 0.5 * d0 - 0.5 * d2;
   // Point3D dirDI_NE = 0.5 * d0 + 0.5 * d2;
   // Point3D dirDI_NW = -dirDI_SE;
   // Point3D dirDI_SW = -dirDI_NE;

   // coords
   Point3D x0  = edge.getCoordinates()[0], x;
   Point3D dx = h * edge.getDirection();

   // forms
   P2Form formS, formN;
   formS.setGeometryMap(faceS->getGeometryMap());

   if (faceN) formN.setGeometryMap(faceN->getGeometryMap());

   real_t tmp = walberla::real_c(0);

   // stencil entries
   std::vector<real_t> vertexToVertexStencil(7);
   std::vector<real_t> edgeToVertexStencil(12);
   std::vector<real_t> vertexToEdgeStencil(12);
   std::vector<real_t> edgeToEdgeStencil(15);

   // loop over all DOFs
   for (const auto& it : hyteg::edgedof::macroedge::Iterator(level, 0))
   {
      uint_t i = it.col();
      std::fill(vertexToVertexStencil.begin(), vertexToVertexStencil.end(), 0.0);
      std::fill(edgeToVertexStencil.begin(), edgeToVertexStencil.end(), 0.0);
      std::fill(edgeToEdgeStencil.begin(), edgeToEdgeStencil.end(), 0.0);
      std::fill(vertexToEdgeStencil.begin(), vertexToEdgeStencil.end(), 0.0);

      ////////// VERTEX //////////
      if (i != 0)
      {
         // assemble stencils
         x = x0 + walberla::real_c(i) * dx;
         /// south face
         //// vertex to vertex
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(formS, {x, x + dir_W, x + dir_S}, P1Elements::P1Elements2D::elementSW, vertexToVertexStencil.data());
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(formS, {x, x + dir_S, x + dir_SE}, P1Elements::P1Elements2D::elementS, vertexToVertexStencil.data());
         vertexdof::variablestencil::assembleLocalStencil<P2Form>(formS, {x, x + dir_SE, x + dir_E}, P1Elements::P1Elements2D::elementSE, vertexToVertexStencil.data());
         //// edge to vertex
         assembleEdgeToVertexStencil<P2Form>(formS, {x, x + dir_W, x + dir_S}, P2Elements::P2Face::elementSW_reord, edgeToVertexStencil.data());
         assembleEdgeToVertexStencil<P2Form>(formS, {x, x + dir_S, x + dir_SE}, P2Elements::P2Face::elementS_reord, edgeToVertexStencil.data());
         assembleEdgeToVertexStencil<P2Form>(formS, {x, x + dir_SE, x + dir_E}, P2Elements::P2Face::elementSE_reord, edgeToVertexStencil.data());

         /// north face
         if (faceN)
         {
            //// vertex to vertex
            vertexdof::variablestencil::assembleLocalStencil<P2Form>(formN, {x, x + dir_E, x + dir_N}, P1Elements::P1Elements2D::elementNE, vertexToVertexStencil.data());
            vertexdof::variablestencil::assembleLocalStencil<P2Form>(formN, {x, x + dir_N, x + dir_NW}, P1Elements::P1Elements2D::elementN, vertexToVertexStencil.data());
            vertexdof::variablestencil::assembleLocalStencil<P2Form>(formN, {x, x + dir_NW, x + dir_W}, P1Elements::P1Elements2D::elementNW, vertexToVertexStencil.data());
            //// edge to vertex
            assembleEdgeToVertexStencil<P2Form>(formN, {x, x + dir_E, x + dir_N}, P2Elements::P2Face::elementNE_reord, edgeToVertexStencil.data());
            assembleEdgeToVertexStencil<P2Form>(formN, {x, x + dir_N, x + dir_NW}, P2Elements::P2Face::elementN_reord, edgeToVertexStencil.data());
            assembleEdgeToVertexStencil<P2Form>(formN, {x, x + dir_NW, x + dir_W}, P2Elements::P2Face::elementNW_reord, edgeToVertexStencil.data());
         }

         // apply GS
         tmp = rhsVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, SD::VERTEX_C)];

         /// on edge vertex dof
         for (const auto& dir : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF)
         {
            tmp -= dstVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, dir)] * vertexToVertexStencil[vertexdof::stencilIndexFromVertex(dir)];
         }

         /// on edge edge dof
         for (const auto& dir : edgedof::macroedge::neighborsOnEdgeFromVertex)
         {
            tmp -= dstEdgeDoF[edgedof::macroedge::indexFromVertex(level, i, dir)] * edgeToVertexStencil[edgedof::stencilIndexFromVertex(dir)];
         }

         /// south face vertex dof
         for (const auto& dir : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF)
         {
            tmp -= dstVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, dir)] * vertexToVertexStencil[vertexdof::stencilIndexFromVertex(dir)];
         }

         /// south face edge dof
         for (const auto& dir : edgedof::macroedge::neighborsOnSouthFaceFromVertex)
         {
            tmp -= dstEdgeDoF[edgedof::macroedge::indexFromVertex(level, i, dir)] * edgeToVertexStencil[edgedof::stencilIndexFromVertex(dir)];
         }

         if (faceN)
         {
            /// north face vertex dof
            for (const auto& dir : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF)
            {
               tmp -= dstVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, dir)] * vertexToVertexStencil[vertexdof::stencilIndexFromVertex(dir)];
            }

            /// north face edge dof
            for (const auto& dir : edgedof::macroedge::neighborsOnNorthFaceFromVertex)
            {
               tmp -= dstEdgeDoF[edgedof::macroedge::indexFromVertex(level, i, dir)] * edgeToVertexStencil[edgedof::stencilIndexFromVertex(dir)];
            }
         }

         dstVertexDoF[vertexdof::macroedge::indexFromVertex(level, i, SD::VERTEX_C)] =
            tmp / vertexToVertexStencil[vertexdof::stencilIndexFromVertex(SD::VERTEX_C)];
      }

      ////////// (HORIZONTAL) EDGE //////////
      // assemble stencils
      x = x0 + walberla::real_c(i + 0.5) * dx;
      /// south face
      //// vertex to edge
      assembleVertexToEdgeStencil<P2Form>(formS, {x + dirHO_W, x + dirHO_E, x + dirHO_SE}, HO_v_S, vertexToEdgeStencil.data());
      //// edge to edge
      assembleEdgeToEdgeStencil<P2Form>(formS, {x + dirHO_W, x + dirHO_E, x + dirHO_SE}, HO_e_S, edgeToEdgeStencil.data());

      /// north face
      if (faceN)
      {
         //// vertex to edge
         assembleVertexToEdgeStencil<P2Form>(formN, {x + dirHO_W, x + dirHO_E, x + dirHO_NW}, HO_v_N, vertexToEdgeStencil.data());
         //// edge to edge
         assembleEdgeToEdgeStencil<P2Form>(formN, {x + dirHO_W, x + dirHO_E, x + dirHO_NW}, HO_e_N, edgeToEdgeStencil.data());
      }

      // apply GS
      tmp = rhsEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge(level, i, SD::EDGE_HO_C)];

      /// on edge
      for (const auto& dir : vertexdof::macroedge::neighborsOnEdgeFromHorizontalEdgeDoF)
      {
         tmp -= dstVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge(level, i, dir)] * vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge(dir)];
      }

      /// on south face
      for (const auto& dir : vertexdof::macroedge::neighborsOnSouthFaceFromHorizontalEdgeDoF)
      {
         tmp -= dstVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge(level, i, dir)] * vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge(dir)];
      }

      for (const auto& dir : edgedof::macroedge::neighborsOnSouthFaceFromHorizontalEdge)
      {
         tmp -= dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge(level, i, dir)] * edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge(dir)];
      }

      /// on north face
      if (faceN)
      {
         for (const auto& dir : vertexdof::macroedge::neighborsOnNorthFaceFromHorizontalEdgeDoF)
         {
            tmp -= dstVertexDoF[vertexdof::macroedge::indexFromHorizontalEdge(level, i, dir)] * vertexToEdgeStencil[vertexdof::stencilIndexFromHorizontalEdge(dir)];
         }

         for (const auto& dir : edgedof::macroedge::neighborsOnNorthFaceFromHorizontalEdge)
         {
            tmp -= dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge(level, i, dir)] * edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge(dir)];
         }
      }

      dstEdgeDoF[edgedof::macroedge::indexFromHorizontalEdge(level, i, SD::EDGE_HO_C)] =
         tmp / edgeToEdgeStencil[edgedof::stencilIndexFromHorizontalEdge(SD::EDGE_HO_C)];;
   }
}

} // namespace macroedge

namespace macrovertex {

template <class P2Form>
inline void applyVariableStencil(uint_t level,
                                 const Vertex& vertex,
                                 const std::shared_ptr<PrimitiveStorage>& storage,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Vertex>& srcVertexDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Vertex>& srcEdgeDoFID,
                                 const PrimitiveDataID<FunctionMemory<real_t>, Vertex>& dstVertexDoFID,
                                 UpdateType update)
{
   typedef stencilDirection SD;

   real_t* srcVertexDoF = vertex.getData(srcVertexDoFID)->getPointer(level);
   real_t* srcEdgeDoF   = vertex.getData(srcEdgeDoFID)->getPointer(level);
   real_t* dstVertexDoF = vertex.getData(dstVertexDoFID)->getPointer(level);

   const uint_t NE = vertex.getNumNeighborEdges();
   const uint_t NF = vertex.getNumNeighborFaces();
   const real_t  h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));

   // assemble stencils
   Point3D x, d0, d2;
   P2Form form;
   std::vector<real_t> vertexToVertexStencil(NE + 1);
   std::vector<real_t> edgeToVertexStencil(NE + NF);
   std::fill(vertexToVertexStencil.begin(), vertexToVertexStencil.end(), 0.0);
   std::fill(edgeToVertexStencil.begin(), edgeToVertexStencil.end(), 0.0);

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

      Point3D vertexToVertex;
      Point3D edgeToVertex;
      form.integrate({{x, x + d0, x + d2}}, vertexToVertex);
      form.integrateEdgeToVertex({{x, x + d0, x + d2}}, edgeToVertex);

      // center contribution
      vertexToVertexStencil[0] += vertexToVertex[0];

      // neighbour edge contribution
      for (uint_t i = 0; i < ne; ++i)
      {
         uint_t edge_idx = vertex.edge_index(adj_edges[i]);
         edgeToVertexStencil[edge_idx] += edgeToVertex[2 - i];
         vertexToVertexStencil[edge_idx + 1] += vertexToVertex[i + 1];
      }

      // opposite edge contribution
      edgeToVertexStencil[face_idx] += edgeToVertex[0];
   }

   // apply Operator
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

template <class P2Form>
inline void smoothGSVariableStencil(uint_t level,
                                    const Vertex& vertex,
                                    const std::shared_ptr<PrimitiveStorage>& storage,
                                    const PrimitiveDataID<FunctionMemory<real_t>, Vertex>& dstVertexDoFID,
                                    const PrimitiveDataID<FunctionMemory<real_t>, Vertex>& dstEdgeDoFID,
                                    const PrimitiveDataID<FunctionMemory<real_t>, Vertex>& rhsVertexDoFID)
{
   typedef stencilDirection SD;

   real_t* dstVertexDoF = vertex.getData(dstVertexDoFID)->getPointer(level);
   real_t* dstEdgeDoF   = vertex.getData(dstEdgeDoFID)->getPointer(level);
   real_t* rhsVertexDoF = vertex.getData(rhsVertexDoFID)->getPointer(level);

   const uint_t NE = vertex.getNumNeighborEdges();
   const uint_t NF = vertex.getNumNeighborFaces();
   const real_t  h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));

   // assemble stencils
   Point3D x, d0, d2;
   P2Form form;
   std::vector<real_t> vertexToVertexStencil(NE + 1);
   std::vector<real_t> edgeToVertexStencil(NE + NF);
   std::fill(vertexToVertexStencil.begin(), vertexToVertexStencil.end(), 0.0);
   std::fill(edgeToVertexStencil.begin(), edgeToVertexStencil.end(), 0.0);

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

      Point3D vertexToVertex;
      Point3D edgeToVertex;
      form.integrate({{x, x + d0, x + d2}}, vertexToVertex);
      form.integrateEdgeToVertex({{x, x + d0, x + d2}}, edgeToVertex);

      // center contribution
      vertexToVertexStencil[0] += vertexToVertex[0];

      // neighbour edge contribution
      for (uint_t i = 0; i < ne; ++i)
      {
         uint_t edge_idx = vertex.edge_index(adj_edges[i]);
         edgeToVertexStencil[edge_idx] += edgeToVertex[2 - i];
         vertexToVertexStencil[edge_idx + 1] += vertexToVertex[i + 1];
      }

      // opposite edge contribution
      edgeToVertexStencil[face_idx] += edgeToVertex[0];
   }

   // apply GS
   dstVertexDoF[0] = rhsVertexDoF[0];

   for (size_t i = 1; i < NE + 1; ++i)
   {
      dstVertexDoF[0] -= vertexToVertexStencil[i] * dstVertexDoF[i];
   }

   for (size_t i = 0; i < NE + NF; ++i)
   {
      dstVertexDoF[0] -= edgeToVertexStencil[i] * dstEdgeDoF[i];
   }

   dstVertexDoF[0] /= vertexToVertexStencil[0];
}

} // namespace macrovertex
} // namespace variablestencil
} // namespace P2
} // namespace hyteg