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

#include "core/debug/all.h"

#include "hyteg/Levelinfo.hpp"
#include "hyteg/Macros.hpp"
#include "hyteg/p1functionspace/P1Elements.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/polynomial/PolynomialEvaluator.hpp"
#include "hyteg/primitives/Face.hpp"

namespace hyteg {
namespace vertexdof {
namespace variablestencil {

template < class P1Form >
inline void assembleLocalStencil( const P1Form&                            form,
                                  const std::array< Point3D, 3 >&          coords,
                                  const std::array< stencilDirection, 3 >& directions,
                                  real_t*                                  opr_data )
{
   Point3D matrixRow;

   form.integrate( coords, matrixRow );

   opr_data[vertexdof::stencilIndexFromVertex( directions[0] )] += matrixRow[0];
   opr_data[vertexdof::stencilIndexFromVertex( directions[1] )] += matrixRow[1];
   opr_data[vertexdof::stencilIndexFromVertex( directions[2] )] += matrixRow[2];
}

// as above but using the new integrateRow()-interface. Old version is kept for legacy purposes, e.g., P2 Operators.
// todo: remove old version once all Operators are renewed
template < class P1Form >
inline void assembleLocalStencil_new( const P1Form&                            form,
                                  const std::array< Point3D, 3 >&          coords,
                                  const std::array< stencilDirection, 3 >& directions,
                                  real_t*                                  opr_data )
{
   Matrixr<1,3> matrixRow;

   form.integrateRow(0, coords, matrixRow );

   opr_data[vertexdof::stencilIndexFromVertex( directions[0] )] += matrixRow(0,0);
   opr_data[vertexdof::stencilIndexFromVertex( directions[1] )] += matrixRow(0,1);
   opr_data[vertexdof::stencilIndexFromVertex( directions[2] )] += matrixRow(0,2);
}

namespace macroface {

template < typename ValueType, class P1Form >
inline void applyVariableStencil(uint_t Level,
                                 const Face &face,
                                 const PrimitiveDataID<FunctionMemory<ValueType>, Face> &srcId,
                                 const PrimitiveDataID<FunctionMemory<ValueType>, Face> &dstId,
                                 UpdateType update)
{
   typedef stencilDirection SD;

   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto src = face.getData( srcId )->getPointer( Level );
   auto dst = face.getData( dstId )->getPointer( Level );

   Point3D x0( face.coords[0] ), x;
   real_t  h = 1.0 / ( walberla::real_c( rowsize - 1 ) );

   Point3D d0 = h * ( face.coords[1] - face.coords[0] );
   Point3D d2 = h * ( face.coords[2] - face.coords[0] );

   P1Form form;
   form.setGeometryMap( face.getGeometryMap() );

   ValueType tmp;

   Point3D dirS  = -1.0 * d2;
   Point3D dirSE = d0 - 1.0 * d2;
   Point3D dirE  = d0;
   Point3D dirW  = -1.0 * d0;
   Point3D dirNW = -1.0 * d0 + d2;
   Point3D dirN  = d2;

   std::vector< real_t > opr_data( 7 );

   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      x = x0;
      x += walberla::real_c( j ) * d2 + d0;

      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         std::fill( opr_data.begin(), opr_data.end(), 0.0 );

         assembleLocalStencil< P1Form >( form, {x, x + dirW, x + dirS}, P1Elements::P1Elements2D::elementSW, opr_data.data() );
         assembleLocalStencil< P1Form >( form, {x, x + dirS, x + dirSE}, P1Elements::P1Elements2D::elementS, opr_data.data() );
         assembleLocalStencil< P1Form >( form, {x, x + dirSE, x + dirE}, P1Elements::P1Elements2D::elementSE, opr_data.data() );
         assembleLocalStencil< P1Form >( form, {x, x + dirE, x + dirN}, P1Elements::P1Elements2D::elementNE, opr_data.data() );
         assembleLocalStencil< P1Form >( form, {x, x + dirN, x + dirNW}, P1Elements::P1Elements2D::elementN, opr_data.data() );
         assembleLocalStencil< P1Form >( form, {x, x + dirNW, x + dirW}, P1Elements::P1Elements2D::elementNW, opr_data.data() );

         //      PointND<real_t, 7> test(opr_data.data());
         //      WALBERLA_LOG_INFO("stencil = " << test);

         if( update == Replace )
         {
            tmp = ValueType( 0 );
         } else
         {
            tmp = dst[vertexdof::macroface::indexFromVertex( Level, i, j, SD::VERTEX_C )];
         }

         tmp += opr_data[vertexdof::stencilIndexFromVertex( SD::VERTEX_C )] *
                src[vertexdof::macroface::indexFromVertex( Level, i, j, SD::VERTEX_C )];
         tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[0] )] *
                src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[0] )];
         tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[1] )] *
                src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[1] )];
         tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[2] )] *
                src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[2] )];
         tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[3] )] *
                src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[3] )];
         tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[4] )] *
                src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[4] )];
         tmp += opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[5] )] *
                src[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[5] )];

         dst[vertexdof::macroface::indexFromVertex( Level, i, j, SD::VERTEX_C )] = tmp;

         x += d0;
      }
      --inner_rowsize;
   }
}

template < typename ValueType, class P1Form >
inline void smoothGSVariableStencil(uint_t Level,
                                    const Face &face,
                                    const PrimitiveDataID<FunctionMemory<ValueType>, Face> &dstId,
                                    const PrimitiveDataID<FunctionMemory<ValueType>, Face> &rhsId)
{
   uint_t rowsize       = levelinfo::num_microvertices_per_edge( Level );
   uint_t inner_rowsize = rowsize;

   auto rhs = face.getData( rhsId )->getPointer( Level );
   auto dst = face.getData( dstId )->getPointer( Level );

   Point3D x0( face.coords[0] ), x;
   real_t  h = 1.0 / ( walberla::real_c( rowsize - 1 ) );

   Point3D d0 = h * ( face.coords[1] - face.coords[0] );
   Point3D d2 = h * ( face.coords[2] - face.coords[0] );

   P1Form form;
   form.setGeometryMap( face.getGeometryMap() );

   ValueType tmp;

   Point3D dirS  = -1.0 * d2;
   Point3D dirSE = d0 - 1.0 * d2;
   Point3D dirE  = d0;
   Point3D dirW  = -1.0 * d0;
   Point3D dirNW = -1.0 * d0 + d2;
   Point3D dirN  = d2;

   std::vector< real_t > opr_data( 7 );


   for( uint_t j = 1; j < rowsize - 2; ++j )
   {
      x = x0;
      x += walberla::real_c( j ) * d2 + d0;

      for( uint_t i = 1; i < inner_rowsize - 2; ++i )
      {
         std::fill( opr_data.begin(), opr_data.end(), 0.0 );

         assembleLocalStencil< P1Form >( form, {x, x + dirW, x + dirS}, P1Elements::P1Elements2D::elementSW, opr_data.data() );
         assembleLocalStencil< P1Form >( form, {x, x + dirS, x + dirSE}, P1Elements::P1Elements2D::elementS, opr_data.data() );
         assembleLocalStencil< P1Form >( form, {x, x + dirSE, x + dirE}, P1Elements::P1Elements2D::elementSE, opr_data.data() );
         assembleLocalStencil< P1Form >( form, {x, x + dirE, x + dirN}, P1Elements::P1Elements2D::elementNE, opr_data.data() );
         assembleLocalStencil< P1Form >( form, {x, x + dirN, x + dirNW}, P1Elements::P1Elements2D::elementN, opr_data.data() );
         assembleLocalStencil< P1Form >( form, {x, x + dirNW, x + dirW}, P1Elements::P1Elements2D::elementNW, opr_data.data() );

         //      PointND<real_t, 7> test(opr_data.data());
         //      WALBERLA_LOG_INFO("stencil = " << test);

         tmp = rhs[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )];

         //for (auto neighbor : neighbors) {
         for( uint_t k = 0; k < vertexdof::macroface::neighborsWithoutCenter.size(); ++k )
         {
            tmp -= opr_data[vertexdof::stencilIndexFromVertex( vertexdof::macroface::neighborsWithoutCenter[k] )] *
                   dst[vertexdof::macroface::indexFromVertex( Level, i, j, vertexdof::macroface::neighborsWithoutCenter[k] )];
         }

         dst[vertexdof::macroface::indexFromVertex( Level, i, j, stencilDirection::VERTEX_C )] =
             tmp / opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];

         x += d0;
      }
      --inner_rowsize;
   }
}

} // namespace macroface

namespace macroedge {

template < typename ValueType, class P1Form >
inline void applyVariableStencil(uint_t Level,
                                 const Edge &edge,
                                 const std::shared_ptr<PrimitiveStorage> &storage,
                                 const PrimitiveDataID<FunctionMemory<ValueType>, Edge> &srcId,
                                 const PrimitiveDataID<FunctionMemory<ValueType>, Edge> &dstId,
                                 UpdateType update)
{
   size_t rowsize = levelinfo::num_microvertices_per_edge( Level );

   auto src = edge.getData( srcId )->getPointer( Level );
   auto dst = edge.getData( dstId )->getPointer( Level );

   ValueType tmp;

   Face*  faceS = storage->getFace( edge.neighborFaces()[0] );
   Face*  faceN = nullptr;
   uint_t s_south = faceS->vertex_index( edge.neighborVertices()[0] );
   uint_t e_south = faceS->vertex_index( edge.neighborVertices()[1] );
   uint_t o_south = faceS->vertex_index( faceS->get_vertex_opposite_to_edge( edge.getID() ) );

   real_t h = 1.0 / ( walberla::real_c( rowsize - 1 ) );

   Point3D dS_se = h * ( faceS->coords[e_south] - faceS->coords[s_south] );
   Point3D dS_so = h * ( faceS->coords[o_south] - faceS->coords[s_south] );
   Point3D dS_oe = h * ( faceS->coords[e_south] - faceS->coords[o_south] );

   Point3D dir_S  = -1.0 * dS_oe;
   Point3D dir_E  = dS_se;
   Point3D dir_SE = dS_so;
   Point3D dir_W  = -1.0 * dS_se;

   Point3D x  = edge.getCoordinates()[0];
   Point3D dx = h * edge.getDirection();
   x += dx;

   uint_t  s_north, e_north, o_north;
   Point3D dir_N;
   Point3D dir_NW;
   P1Form form;

   if( edge.getNumNeighborFaces() == 2 )
   {
      faceN   = storage->getFace( edge.neighborFaces()[1] );
      s_north = faceN->vertex_index( edge.neighborVertices()[0] );
      e_north = faceN->vertex_index( edge.neighborVertices()[1] );
      o_north = faceN->vertex_index( faceN->get_vertex_opposite_to_edge( edge.getID() ) );

      Point3D dN_so = h * ( faceN->coords[o_north] - faceN->coords[s_north] );
      Point3D dN_oe = h * ( faceN->coords[e_north] - faceN->coords[o_north] );

      dir_N  = dN_so;
      dir_NW = -1.0 * dN_oe;
   }

   std::vector< real_t > opr_data( 7 );

   for( size_t i = 1; i < rowsize - 1; ++i )
   {
      std::fill( opr_data.begin(), opr_data.end(), 0.0 );

      // assemble south
      form.setGeometryMap( faceS->getGeometryMap() );
      assembleLocalStencil< P1Form >( form, {x, x + dir_W, x + dir_S}, P1Elements::P1Elements2D::elementSW, opr_data.data() );
      assembleLocalStencil< P1Form >( form, {x, x + dir_S, x + dir_SE}, P1Elements::P1Elements2D::elementS, opr_data.data() );
      assembleLocalStencil< P1Form >( form, {x, x + dir_SE, x + dir_E}, P1Elements::P1Elements2D::elementSE, opr_data.data() );

      if( edge.getNumNeighborFaces() == 2 )
      {
         form.setGeometryMap( faceN->getGeometryMap() );
         assembleLocalStencil< P1Form >( form, {x, x + dir_E, x + dir_N}, P1Elements::P1Elements2D::elementNE, opr_data.data() );
         assembleLocalStencil< P1Form >( form, {x, x + dir_N, x + dir_NW}, P1Elements::P1Elements2D::elementN, opr_data.data() );
         assembleLocalStencil< P1Form >( form, {x, x + dir_NW, x + dir_W}, P1Elements::P1Elements2D::elementNW, opr_data.data() );
      }

      tmp = opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] *
            src[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )];

      // neighbors on edge
      for( const auto& neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
      {
         tmp += opr_data[vertexdof::stencilIndexFromVertex( neighbor )] *
                src[vertexdof::macroedge::indexFromVertex( Level, i, neighbor )];
      }

      for( const auto& neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
      {
         tmp += opr_data[vertexdof::stencilIndexFromVertex( neighbor )] *
                src[vertexdof::macroedge::indexFromVertex( Level, i, neighbor )];
      }

      if( edge.getNumNeighborFaces() == 2 )
      {
         for( const auto& neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
         {
            tmp += opr_data[vertexdof::stencilIndexFromVertex( neighbor )] *
                   src[vertexdof::macroedge::indexFromVertex( Level, i, neighbor )];
         }
      }

      if( update == Replace )
      {
         dst[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] = tmp;
      } else if( update == Add )
      {
         dst[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] += tmp;
      }

      x += dx;
   }
}

template < typename ValueType, class P1Form >
inline void smoothGSVariableStencil(uint_t Level,
                                    const Edge &edge,
                                    const std::shared_ptr<PrimitiveStorage> &storage,
                                    const PrimitiveDataID<FunctionMemory<ValueType>, Edge> &dstId,
                                    const PrimitiveDataID<FunctionMemory<ValueType>, Edge> &rhsId)
{
   size_t rowsize = levelinfo::num_microvertices_per_edge( Level );

   auto rhs = edge.getData( rhsId )->getPointer( Level );
   auto dst = edge.getData( dstId )->getPointer( Level );

   Face*  faceS = storage->getFace( edge.neighborFaces()[0] );
   Face*  faceN = nullptr;
   uint_t s_south = faceS->vertex_index( edge.neighborVertices()[0] );
   uint_t e_south = faceS->vertex_index( edge.neighborVertices()[1] );
   uint_t o_south = faceS->vertex_index( faceS->get_vertex_opposite_to_edge( edge.getID() ) );

   real_t h = 1.0 / ( walberla::real_c( rowsize - 1 ) );

   Point3D dS_se = h * ( faceS->coords[e_south] - faceS->coords[s_south] );
   Point3D dS_so = h * ( faceS->coords[o_south] - faceS->coords[s_south] );
   Point3D dS_oe = h * ( faceS->coords[e_south] - faceS->coords[o_south] );

   Point3D dir_S  = -1.0 * dS_oe;
   Point3D dir_E  = dS_se;
   Point3D dir_SE = dS_so;
   Point3D dir_W  = -1.0 * dS_se;

   Point3D x  = edge.getCoordinates()[0];
   Point3D dx = h * edge.getDirection();
   x += dx;

   uint_t  s_north, e_north, o_north;
   Point3D dir_N;
   Point3D dir_NW;
   P1Form form;

   if( edge.getNumNeighborFaces() == 2 )
   {
      faceN   = storage->getFace( edge.neighborFaces()[1] );
      s_north = faceN->vertex_index( edge.neighborVertices()[0] );
      e_north = faceN->vertex_index( edge.neighborVertices()[1] );
      o_north = faceN->vertex_index( faceN->get_vertex_opposite_to_edge( edge.getID() ) );

      Point3D dN_so = h * ( faceN->coords[o_north] - faceN->coords[s_north] );
      Point3D dN_oe = h * ( faceN->coords[e_north] - faceN->coords[o_north] );

      dir_N  = dN_so;
      dir_NW = -1.0 * dN_oe;
   }

   std::vector< real_t > opr_data( 7 );

   for( size_t i = 1; i < rowsize - 1; ++i )
   {
      std::fill( opr_data.begin(), opr_data.end(), 0.0 );

      // assemble south
      form.setGeometryMap( faceS->getGeometryMap() );
      assembleLocalStencil< P1Form >( form, {x, x + dir_W, x + dir_S}, P1Elements::P1Elements2D::elementSW, opr_data.data() );
      assembleLocalStencil< P1Form >( form, {x, x + dir_S, x + dir_SE}, P1Elements::P1Elements2D::elementS, opr_data.data() );
      assembleLocalStencil< P1Form >( form, {x, x + dir_SE, x + dir_E}, P1Elements::P1Elements2D::elementSE, opr_data.data() );

      if( edge.getNumNeighborFaces() == 2 )
      {
         form.setGeometryMap( faceN->getGeometryMap() );
         assembleLocalStencil< P1Form >( form, {x, x + dir_E, x + dir_N}, P1Elements::P1Elements2D::elementNE, opr_data.data() );
         assembleLocalStencil< P1Form >( form, {x, x + dir_N, x + dir_NW}, P1Elements::P1Elements2D::elementN, opr_data.data() );
         assembleLocalStencil< P1Form >( form, {x, x + dir_NW, x + dir_W}, P1Elements::P1Elements2D::elementNW, opr_data.data() );
      }

      dst[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] =
          rhs[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )];

      for( const auto& neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
      {
         dst[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] -=
             opr_data[vertexdof::stencilIndexFromVertex( neighbor )] *
             dst[vertexdof::macroedge::indexFromVertex( Level, i, neighbor )];
      }

      for( const auto& neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
      {
         dst[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] -=
             opr_data[vertexdof::stencilIndexFromVertex( neighbor )] *
             dst[vertexdof::macroedge::indexFromVertex( Level, i, neighbor )];
      }

      if( edge.getNumNeighborFaces() == 2 )
      {
         for( const auto& neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
         {
            dst[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] -=
                opr_data[vertexdof::stencilIndexFromVertex( neighbor )] *
                dst[vertexdof::macroedge::indexFromVertex( Level, i, neighbor )];
         }
      }

      dst[vertexdof::macroedge::indexFromVertex( Level, i, stencilDirection::VERTEX_C )] /=
          opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];

      x += dx;
   }
}
} // namespace macroedge

namespace macrovertex {

template < typename ValueType, class P1Form >
inline void applyVariableStencil(uint_t level,
                                 const Vertex &vertex,
                                 const std::shared_ptr<PrimitiveStorage> &storage,
                                 const PrimitiveDataID<FunctionMemory<ValueType>, Vertex> &srcId,
                                 const PrimitiveDataID<FunctionMemory<ValueType>, Vertex> &dstId,
                                 UpdateType update)
{
   auto src = vertex.getData( srcId )->getPointer( level );
   auto dst = vertex.getData( dstId )->getPointer( level );

   uint_t rowsize = levelinfo::num_microvertices_per_edge( level );

   std::vector< real_t > opr_data( 1 + vertex.getNumNeighborEdges() );
   std::fill( opr_data.begin(), opr_data.end(), 0.0 );
   Point3D x;
   Point3D d0;
   Point3D d2;
   P1Form form;

   real_t h = 1.0 / ( walberla::real_c( rowsize - 1 ) );

   uint_t neighborId = 0;
   for( auto& faceId : vertex.neighborFaces() )
   {
      Face* face       = storage->getFace( faceId );
      form.setGeometryMap( face->getGeometryMap() );

      uint_t                     v_i       = face->vertex_index( vertex.getID() );
      std::vector< PrimitiveID > adj_edges = face->adjacent_edges( vertex.getID() );

      x = face->coords[v_i];
      d0 =
          ( face->coords[face->vertex_index( storage->getEdge( adj_edges[0] )->get_opposite_vertex( vertex.getID() ) )] - x ) * h;
      d2 =
          ( face->coords[face->vertex_index( storage->getEdge( adj_edges[1] )->get_opposite_vertex( vertex.getID() ) )] - x ) * h;

      Point3D matrixRow;
      form.integrate( {{x, x + d0, x + d2}}, matrixRow );

      uint_t i = 1;
      // iterate over adjacent edges
      for( auto& edgeId : adj_edges )
      {
         uint_t      edge_idx = vertex.edge_index( edgeId ) + 1;
         Edge*       edge     = storage->getEdge( edgeId );
         PrimitiveID vertex_j = edge->get_opposite_vertex( vertex.getID() );

         opr_data[edge_idx] += matrixRow[i];
         i += 1;
      }

      // add contribution of center vertex
      opr_data[0] += matrixRow[0];

      ++neighborId;
   }

   if( update == Replace )
   {
      dst[0] = opr_data[0] * src[0];
   } else if( update == Add )
   {
      dst[0] += opr_data[0] * src[0];
   }

   for( size_t i = 0; i < vertex.getNumNeighborEdges(); ++i )
   {
      dst[0] += opr_data[i + 1] * src[i + 1];
   }
}

template < typename ValueType, class P1Form >
inline void smoothGSVariableStencil(uint_t level,
                                    Vertex &vertex,
                                    const std::shared_ptr<PrimitiveStorage> &storage,
                                    const PrimitiveDataID<FunctionMemory<ValueType>, Vertex> &dstId,
                                    const PrimitiveDataID<FunctionMemory<ValueType>, Vertex> &rhsId)
{
   auto rhs = vertex.getData( rhsId )->getPointer( level );
   auto dst = vertex.getData( dstId )->getPointer( level );

   uint_t rowsize = levelinfo::num_microvertices_per_edge( level );

   std::vector< real_t > opr_data( 1 + vertex.getNumNeighborEdges() );
   std::fill( opr_data.begin(), opr_data.end(), 0.0 );
   Point3D x;
   Point3D d0;
   Point3D d2;

   real_t h = 1.0 / ( walberla::real_c( rowsize - 1 ) );
   P1Form form;

   uint_t neighborId = 0;
   for( auto& faceId : vertex.neighborFaces() )
   {
      Face* face       = storage->getFace( faceId );
      form.setGeometryMap( face->getGeometryMap() );

      uint_t                     v_i       = face->vertex_index( vertex.getID() );
      std::vector< PrimitiveID > adj_edges = face->adjacent_edges( vertex.getID() );

      x = face->coords[v_i];
      d0 =
          ( face->coords[face->vertex_index( storage->getEdge( adj_edges[0] )->get_opposite_vertex( vertex.getID() ) )] - x ) * h;
      d2 =
          ( face->coords[face->vertex_index( storage->getEdge( adj_edges[1] )->get_opposite_vertex( vertex.getID() ) )] - x ) * h;

      Point3D matrixRow;
      form.integrate( {{x, x + d0, x + d2}}, matrixRow );

      uint_t i = 1;
      // iterate over adjacent edges
      for( auto& edgeId : adj_edges )
      {
         uint_t      edge_idx = vertex.edge_index( edgeId ) + 1;
         Edge*       edge     = storage->getEdge( edgeId );
         PrimitiveID vertex_j = edge->get_opposite_vertex( vertex.getID() );

         opr_data[edge_idx] += matrixRow[i];
         i += 1;
      }

      // add contribution of center vertex
      opr_data[0] += matrixRow[0];

      ++neighborId;
   }

   dst[0] = rhs[0];

   for( size_t i = 0; i < vertex.getNumNeighborEdges(); ++i )
   {
      dst[0] -= opr_data[i + 1] * dst[i + 1];
   }

   dst[0] /= opr_data[0];
}

} // namespace macrovertex

} // namespace variablestencil
} // namespace vertexdof
} // namespace hyteg
