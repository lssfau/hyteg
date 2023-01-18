/*
 * Copyright (c) 2017-2023 Daniel Drzisga, Dominik Thoennes, Nils Kohl,
 * Marcus Mohr.
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

#include <array>

#include "hyteg/dgfunctionspace_old/DGFunction.hpp"
#include "hyteg/dgfunctionspace_old/DGFunctionMacroEdge.hpp"
#include "hyteg/dgfunctionspace_old/DGFunctionMacroFace.hpp"
#include "hyteg/dgfunctionspace_old/DGFunctionMacroVertex.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

using walberla::real_c;

class DG0P1UpwindOperator : public Operator< DGFunction_old< real_t >, DGFunction_old< real_t > >
{
 public:
   DG0P1UpwindOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                        const P1VectorFunction< real_t >&          velocity,
                        uint_t                                     minLevel,
                        uint_t                                     maxLevel )
   : Operator( storage, minLevel, maxLevel )
   , velocity_( velocity )
   {
      WALBERLA_ASSERT( velocity.getDimension() == 2 );
   }

   ~DG0P1UpwindOperator() override = default;

   // #define DG0_P1_UPWIND_OLD
   void apply( const DGFunction_old< real_t >& src,
               const DGFunction_old< real_t >& dst,
               uint_t                          level,
               DoFType                         flag,
               UpdateType                      updateType = Replace ) const
   {
      // start pulling edge halos
      src.startCommunication< Face, Edge >( level );

      // end pulling edge halos
      src.endCommunication< Face, Edge >( level );

      // start pulling vertex halos
      src.startCommunication< Edge, Vertex >( level );

      // end pulling vertex halos
      src.endCommunication< Edge, Vertex >( level );

      velocity_[0].startCommunication< Edge, Vertex >( level );
      velocity_[1].startCommunication< Edge, Vertex >( level );

      velocity_[0].startCommunication< Face, Edge >( level );
      velocity_[1].startCommunication< Face, Edge >( level );

      velocity_[0].endCommunication< Edge, Vertex >( level );
      velocity_[1].endCommunication< Edge, Vertex >( level );

      for ( auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if ( testFlag( vertexBC, flag ) )
         {
            treatMacroVertex( level, vertex, storage_, src.getVertexDataID(), dst.getVertexDataID(), updateType );
         }
      }

      dst.startCommunication< Vertex, Edge >( level );

      velocity_[0].endCommunication< Face, Edge >( level );
      velocity_[1].endCommunication< Face, Edge >( level );

      for ( auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            // dgfunction::macroedge::upwind< real_t >( level,
            //                                          edge,
            //                                          storage_,
            //                                          src.getEdgeDataID(),
            //                                          dst.getEdgeDataID(),
            //                                          std::array< PrimitiveDataID< FunctionMemory< real_t >, Edge >, 2 >{
            //                                              { velocity_[0].getEdgeDataID(), velocity_[1].getEdgeDataID() } },
            //                                          updateType );
            treatMacroEdge( level, edge, storage_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
         }
      }

      dst.endCommunication< Vertex, Edge >( level );

      dst.startCommunication< Edge, Face >( level );

      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            treatMacroFace( level, face, storage_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
         }
      }

      dst.endCommunication< Edge, Face >( level );
   }

 private:
   P1VectorFunction< real_t > velocity_;

   void treatMacroFace( const uint_t&                                            level,
                        Face&                                                    face,
                        const std::shared_ptr< PrimitiveStorage >&               storage,
                        const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcId,
                        const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                        UpdateType                                               updateType ) const
   {
      using namespace vertexdof::macroface;

      size_t rowsize       = levelinfo::num_microvertices_per_edge( level );
      size_t inner_rowsize = rowsize;

      // get memories
      auto src = face.getData( srcId )->getPointer( level );
      auto dst = face.getData( dstId )->getPointer( level );
      auto u   = face.getData( velocity_[0].getFaceDataID() )->getPointer( level );
      auto v   = face.getData( velocity_[1].getFaceDataID() )->getPointer( level );

      // get edge directions
      auto d0 = ( face.getCoordinates()[1] - face.getCoordinates()[0] ) / real_c( rowsize - 1 );
      auto d1 = ( face.getCoordinates()[2] - face.getCoordinates()[1] ) / real_c( rowsize - 1 );
      auto d2 = ( face.getCoordinates()[0] - face.getCoordinates()[2] ) / real_c( rowsize - 1 );

      // compute edge lengths
      real_t d0Length = d0.norm();
      real_t d1Length = d1.norm();
      real_t d2Length = d2.norm();

      // compute rightward pointing normals
      auto n_0 = d0.normal2D() / d0Length;
      auto n_1 = d1.normal2D() / d1Length;
      auto n_2 = d2.normal2D() / d2Length;

      real_t faceOrientation =
          math::faceOrientation2D( face.getCoordinates()[0], face.getCoordinates()[1], face.getCoordinates()[2] );

      // correct normals if all normals point in wrong direction, i.e. inwards
      n_0 *= faceOrientation;
      n_1 *= faceOrientation;
      n_2 *= faceOrientation;

      real_t faceArea    = std::pow( real_c( 4 ), -real_c( level ) ) * face.getArea();
      real_t faceAreaInv = real_c( 1 ) / faceArea;

      real_t tmp;

      Point2D u_0, u_1, u_2;
      real_t  un_0, un_1, un_2;
      real_t  c_up_0, c_up_1, c_up_2;

      for ( size_t j = 1; j < rowsize - 2; ++j )
      {
         for ( size_t i = 1; i < inner_rowsize - 3; ++i )
         {
            // evaluate velocities
            u_0[0] = real_c( 0.5 ) * ( u[indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] +
                                       u[indexFromVertex( level, i + 1, j, stencilDirection::VERTEX_C )] );
            u_0[1] = real_c( 0.5 ) * ( v[indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] +
                                       v[indexFromVertex( level, i + 1, j, stencilDirection::VERTEX_C )] );

            u_1[0] = real_c( 0.5 ) * ( u[indexFromVertex( level, i + 1, j, stencilDirection::VERTEX_C )] +
                                       u[indexFromVertex( level, i, j + 1, stencilDirection::VERTEX_C )] );
            u_1[1] = real_c( 0.5 ) * ( v[indexFromVertex( level, i + 1, j, stencilDirection::VERTEX_C )] +
                                       v[indexFromVertex( level, i, j + 1, stencilDirection::VERTEX_C )] );

            u_2[0] = real_c( 0.5 ) * ( u[indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] +
                                       u[indexFromVertex( level, i, j + 1, stencilDirection::VERTEX_C )] );
            u_2[1] = real_c( 0.5 ) * ( v[indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] +
                                       v[indexFromVertex( level, i, j + 1, stencilDirection::VERTEX_C )] );

            // compute magnitude of outward normal component and scale be edge length
            un_0 = d0Length * u_0.dot( n_0 );
            un_1 = d1Length * u_1.dot( n_1 );
            un_2 = d2Length * u_2.dot( n_2 );

            if ( un_0 >= 0 )
            {
               c_up_0 = src[facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_GRAY_C )];
            }
            else
            {
               c_up_0 = src[facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_BLUE_S )];
            }

            if ( un_1 >= 0 )
            {
               c_up_1 = src[facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_GRAY_C )];
            }
            else
            {
               c_up_1 = src[facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_BLUE_E )];
            }

            if ( un_2 >= 0 )
            {
               c_up_2 = src[facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_GRAY_C )];
            }
            else
            {
               c_up_2 = src[facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_BLUE_W )];
            }

            tmp = un_0 * c_up_0 + un_1 * c_up_1 + un_2 * c_up_2;
            tmp *= faceAreaInv;

            if ( updateType == Replace )
            {
               dst[facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_GRAY_C )] = tmp;
            }
            else if ( updateType == Add )
            {
               dst[facedof::macroface::indexFaceFromGrayFace( level, i, j, stencilDirection::CELL_GRAY_C )] += tmp;
            }
         }
         --inner_rowsize;
      }

      inner_rowsize = rowsize;

      // flip normals
      n_0 *= -1.0;
      n_1 *= -1.0;
      n_2 *= -1.0;

      for ( size_t j = 0; j < rowsize - 2; ++j )
      {
         for ( size_t i = 0; i < inner_rowsize - 2; ++i )
         {
            // evalate velocities
            u_0[0] = real_c( 0.5 ) * ( u[indexFromVertex( level, i, j + 1, stencilDirection::VERTEX_C )] +
                                       u[indexFromVertex( level, i + 1, j + 1, stencilDirection::VERTEX_C )] );
            u_0[1] = real_c( 0.5 ) * ( v[indexFromVertex( level, i, j + 1, stencilDirection::VERTEX_C )] +
                                       v[indexFromVertex( level, i + 1, j + 1, stencilDirection::VERTEX_C )] );

            u_1[0] = real_c( 0.5 ) * ( u[indexFromVertex( level, i, j + 1, stencilDirection::VERTEX_C )] +
                                       u[indexFromVertex( level, i + 1, j, stencilDirection::VERTEX_C )] );
            u_1[1] = real_c( 0.5 ) * ( v[indexFromVertex( level, i, j + 1, stencilDirection::VERTEX_C )] +
                                       v[indexFromVertex( level, i + 1, j, stencilDirection::VERTEX_C )] );

            u_2[0] = real_c( 0.5 ) * ( u[indexFromVertex( level, i + 1, j, stencilDirection::VERTEX_C )] +
                                       u[indexFromVertex( level, i + 1, j + 1, stencilDirection::VERTEX_C )] );
            u_2[1] = real_c( 0.5 ) * ( v[indexFromVertex( level, i + 1, j, stencilDirection::VERTEX_C )] +
                                       v[indexFromVertex( level, i + 1, j + 1, stencilDirection::VERTEX_C )] );

            un_0 = d0Length * u_0.dot( n_0 );
            un_1 = d1Length * u_1.dot( n_1 );
            un_2 = d2Length * u_2.dot( n_2 );

            if ( un_0 >= 0 )
            {
               c_up_0 = src[facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_BLUE_C )];
            }
            else
            {
               c_up_0 = src[facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_GRAY_N )];
            }

            if ( un_1 >= 0 )
            {
               c_up_1 = src[facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_BLUE_C )];
            }
            else
            {
               c_up_1 = src[facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_GRAY_W )];
            }

            if ( un_2 >= 0 )
            {
               c_up_2 = src[facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_BLUE_C )];
            }
            else
            {
               c_up_2 = src[facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_GRAY_E )];
            }

            tmp = un_0 * c_up_0 + un_1 * c_up_1 + un_2 * c_up_2;
            tmp *= faceAreaInv;

            if ( updateType == Replace )
            {
               dst[facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_BLUE_C )] = tmp;
            }
            else if ( updateType == Add )
            {
               dst[facedof::macroface::indexFaceFromBlueFace( level, i, j, stencilDirection::CELL_BLUE_C )] += tmp;
            }
         }
         --inner_rowsize;
      }
   };

   void treatMacroEdge( const uint_t&                                            level,
                        Edge&                                                    edge,
                        const std::shared_ptr< PrimitiveStorage >&               storage,
                        const PrimitiveDataID< FunctionMemory< real_t >, Edge >& srcId,
                        const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstId,
                        UpdateType                                               updateType ) const
   {
      auto src = edge.getData( srcId )->getPointer( level );
      auto dst = edge.getData( dstId )->getPointer( level );
      auto u   = edge.getData( velocity_[0].getEdgeDataID() )->getPointer( level );
      auto v   = edge.getData( velocity_[1].getEdgeDataID() )->getPointer( level );

      size_t  rowsize = levelinfo::num_microvertices_per_edge( level );
      real_t  tmp;
      Point2D u_0, u_1, u_2;
      real_t  un_0, un_1, un_2;
      real_t  c_up_0, c_up_1, c_up_2;

      // first face (south)
      {
         Face*  face        = storage->getFace( edge.neighborFaces()[0] );
         real_t faceArea    = std::pow( real_c( 4 ), -real_c( level ) ) * face->getArea();
         real_t faceAreaInv = real_c( 1 ) / faceArea;

         auto oppositeVertex = face->get_vertex_opposite_to_edge( edge.getID() );

         uint_t v0 = face->vertex_index( edge.getVertexID0() );
         uint_t v1 = face->vertex_index( edge.getVertexID1() );
         uint_t v2 = face->vertex_index( oppositeVertex );

         // get edge directions
         auto d0 = ( face->getCoordinates()[v1] - face->getCoordinates()[v0] ) / real_c( rowsize - 1 );
         auto d1 = ( face->getCoordinates()[v2] - face->getCoordinates()[v1] ) / real_c( rowsize - 1 );
         auto d2 = ( face->getCoordinates()[v0] - face->getCoordinates()[v2] ) / real_c( rowsize - 1 );

         // compute edge lengths
         real_t d0Length = d0.norm();
         real_t d1Length = d1.norm();
         real_t d2Length = d2.norm();

         // compute normals
         auto n_0 = d0.normal2D() / d0Length;
         auto n_1 = d1.normal2D() / d1Length;
         auto n_2 = d2.normal2D() / d2Length;

         real_t faceOrientation =
             math::faceOrientation2D( face->getCoordinates()[v0], face->getCoordinates()[v1], face->getCoordinates()[v2] );
         n_0 *= faceOrientation;
         n_1 *= faceOrientation;
         n_2 *= faceOrientation;

         for ( uint_t i = 1; i < rowsize - 2; ++i )
         {
            u_0[0] = real_c( 0.5 ) * ( u[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] +
                                       u[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_E )] );
            u_0[1] = real_c( 0.5 ) * ( v[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] +
                                       v[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_E )] );

            u_1[0] = real_c( 0.5 ) * ( u[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_E )] +
                                       u[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_SE )] );
            u_1[1] = real_c( 0.5 ) * ( v[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_E )] +
                                       v[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_SE )] );

            u_2[0] = real_c( 0.5 ) * ( u[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] +
                                       u[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_SE )] );
            u_2[1] = real_c( 0.5 ) * ( v[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] +
                                       v[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_SE )] );

            // CENTER <-- CELL_GRAY_SE (i)
            // NORTH  <-- CELL_GRAY_NE (i)
            // WEST   <-- CELL_BLUE_SE (i)
            // EAST   <-- CELL_BLUE_SE (i+1)

            un_0 = d0Length * u_0.dot( n_0 );
            un_1 = d1Length * u_1.dot( n_1 );
            un_2 = d2Length * u_2.dot( n_2 );

            if ( un_0 >= 0 )
            {
               c_up_0 = src[facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_SE )];
            }
            else
            {
               c_up_0 = src[facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_NE )];
            }

            if ( un_1 >= 0 )
            {
               c_up_1 = src[facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_SE )];
            }
            else
            {
               c_up_1 = src[facedof::macroedge::indexFaceFromVertex( level, i + 1, stencilDirection::CELL_BLUE_SE )];
            }

            if ( un_2 >= 0 )
            {
               c_up_2 = src[facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_SE )];
            }
            else
            {
               c_up_2 = src[facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_BLUE_SE )];
            }

            tmp = un_0 * c_up_0 + un_1 * c_up_1 + un_2 * c_up_2;
            tmp *= faceAreaInv;

            if ( updateType == Replace )
            {
               dst[facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_SE )] = tmp;
            }
            else if ( updateType == Add )
            {
               dst[facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_SE )] += tmp;
            }
         }
      }

      // second face (north)
      if ( edge.getNumNeighborFaces() == 2 )
      {
         Face*  face        = storage->getFace( edge.neighborFaces()[1] );
         real_t faceArea    = std::pow( 4.0, -real_c( level ) ) * face->getArea();
         real_t faceAreaInv = 1.0 / faceArea;

         auto oppositeVertex = face->get_vertex_opposite_to_edge( edge.getID() );

         uint_t v0 = face->vertex_index( edge.getVertexID0() );
         uint_t v1 = face->vertex_index( edge.getVertexID1() );
         uint_t v2 = face->vertex_index( oppositeVertex );

         // get edge directions
         auto d0 = ( face->getCoordinates()[v1] - face->getCoordinates()[v0] ) / real_c( rowsize - 1 );
         auto d1 = ( face->getCoordinates()[v2] - face->getCoordinates()[v1] ) / real_c( rowsize - 1 );
         auto d2 = ( face->getCoordinates()[v0] - face->getCoordinates()[v2] ) / real_c( rowsize - 1 );

         // compute edge lengths
         real_t d0Length = d0.norm();
         real_t d1Length = d1.norm();
         real_t d2Length = d2.norm();

         // compute normals
         auto n_0 = d0.normal2D() / d0Length;
         auto n_1 = d1.normal2D() / d1Length;
         auto n_2 = d2.normal2D() / d2Length;

         real_t faceOrientation =
             math::faceOrientation2D( face->getCoordinates()[v0], face->getCoordinates()[v1], face->getCoordinates()[v2] );
         n_0 *= faceOrientation;
         n_1 *= faceOrientation;
         n_2 *= faceOrientation;

         for ( uint_t i = 1; i < rowsize - 2; ++i )
         {
            u_0[0] = real_c( 0.5 ) * ( u[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] +
                                       u[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_E )] );
            u_0[1] = real_c( 0.5 ) * ( v[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] +
                                       v[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_E )] );

            u_1[0] = real_c( 0.5 ) * ( u[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_E )] +
                                       u[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_N )] );
            u_1[1] = real_c( 0.5 ) * ( v[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_E )] +
                                       v[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_N )] );

            u_2[0] = real_c( 0.5 ) * ( u[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] +
                                       u[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_N )] );
            u_2[1] = real_c( 0.5 ) * ( v[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] +
                                       v[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_N )] );

            // CENTER <-- CELL_GRAY_NE (i)
            // SOUTH  <-- CELL_GRAY_SE (i)
            // WEST   <-- CELL_BLUE_NW (i+1)
            // EAST   <-- CELL_BLUE_NW (i)

            un_0 = d0Length * u_0.dot( n_0 );
            un_1 = d1Length * u_1.dot( n_1 );
            un_2 = d2Length * u_2.dot( n_2 );

            if ( un_0 >= 0 )
            {
               c_up_0 = src[facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_NE )];
            }
            else
            {
               c_up_0 = src[facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_SE )];
            }

            if ( un_1 >= 0 )
            {
               c_up_1 = src[facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_NE )];
            }
            else
            {
               c_up_1 = src[facedof::macroedge::indexFaceFromVertex( level, i + 1, stencilDirection::CELL_BLUE_NW )];
            }

            if ( un_2 >= 0 )
            {
               c_up_2 = src[facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_NE )];
            }
            else
            {
               c_up_2 = src[facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_BLUE_NW )];
            }

            tmp = un_0 * c_up_0 + un_1 * c_up_1 + un_2 * c_up_2;
            tmp *= faceAreaInv;

            if ( updateType == Replace )
            {
               dst[facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_NE )] = tmp;
            }
            else if ( updateType == Add )
            {
               dst[facedof::macroedge::indexFaceFromVertex( level, i, stencilDirection::CELL_GRAY_NE )] += tmp;
            }
         }
      }
   }

   void treatMacroVertex( const uint_t&                                              level,
                          Vertex&                                                    vertex,
                          const std::shared_ptr< PrimitiveStorage >&                 storage,
                          const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& srcId,
                          const PrimitiveDataID< FunctionMemory< real_t >, Vertex >& dstId,
                          UpdateType                                                 updateType ) const
   {
      auto src = vertex.getData( srcId )->getPointer( level );
      auto dst = vertex.getData( dstId )->getPointer( level );
      auto u   = vertex.getData( velocity_[0].getVertexDataID() )->getPointer( level );
      auto v   = vertex.getData( velocity_[1].getVertexDataID() )->getPointer( level );

      size_t  rowsize = levelinfo::num_microvertices_per_edge( level );
      real_t  tmp;
      Point2D u_0, u_1, u_2;
      real_t  un_0, un_1, un_2;
      real_t  c_up_0, c_up_1, c_up_2;

      for ( auto faceIt : vertex.neighborFaces() )
      {
         Face* face = storage->getFace( faceIt );

         real_t faceArea    = std::pow( real_c( 4 ), -real_c( level ) ) * face->getArea();
         real_t faceAreaInv = real_c( 1 ) / faceArea;

         uint_t localFaceId = vertex.face_index( face->getID() );

         uint_t faceMemoryIndex = 2 * localFaceId;
         uint_t blueMemoryIndex = faceMemoryIndex + 1;

         std::vector< PrimitiveID > adjEdgeIds = face->adjacent_edges( vertex.getID() );
         std::vector< Edge* >       adjEdges;
         adjEdges.push_back( storage->getEdge( adjEdgeIds[0] ) );
         adjEdges.push_back( storage->getEdge( adjEdgeIds[1] ) );

         uint_t v0 = face->vertex_index( vertex.getID() );
         uint_t v1 = face->vertex_index( adjEdges[0]->get_opposite_vertex( vertex.getID() ) );
         uint_t v2 = face->vertex_index( adjEdges[1]->get_opposite_vertex( vertex.getID() ) );

         uint_t p1EdgeId0 = vertex.edge_index( adjEdgeIds[0] ) + 1;
         uint_t p1EdgeId1 = vertex.edge_index( adjEdgeIds[1] ) + 1;

         // compute edge directions
         auto d0 = ( face->getCoordinates()[v1] - face->getCoordinates()[v0] ) / real_c( rowsize - 1 );
         auto d1 = ( face->getCoordinates()[v0] - face->getCoordinates()[v2] ) / real_c( rowsize - 1 );
         auto d2 = ( face->getCoordinates()[v2] - face->getCoordinates()[v1] ) / real_c( rowsize - 1 );

         // compute edge lengths
         real_t d0Length = d0.norm();
         real_t d1Length = d1.norm();
         real_t d2Length = d2.norm();

         // compute normals
         auto n_0 = d0.normal2D() / d0Length;
         auto n_1 = d1.normal2D() / d1Length;
         auto n_2 = d2.normal2D() / d2Length;

         real_t faceOrientation =
             math::faceOrientation2D( face->getCoordinates()[v0], face->getCoordinates()[v1], face->getCoordinates()[v2] );
         n_0 *= faceOrientation;
         n_1 *= faceOrientation;
         n_2 *= faceOrientation;

         u_0[0] = real_c( real_c( 0.5 ) ) * ( u[0] + u[p1EdgeId0] );
         u_0[1] = real_c( real_c( 0.5 ) ) * ( v[0] + v[p1EdgeId0] );

         u_1[0] = real_c( real_c( 0.5 ) ) * ( u[0] + u[p1EdgeId1] );
         u_1[1] = real_c( real_c( 0.5 ) ) * ( v[0] + v[p1EdgeId1] );

         u_2[0] = real_c( real_c( 0.5 ) ) * ( u[p1EdgeId0] + u[p1EdgeId1] );
         u_2[1] = real_c( real_c( 0.5 ) ) * ( v[p1EdgeId0] + v[p1EdgeId1] );

         un_0 = d0Length * u_0.dot( n_0 );
         un_1 = d1Length * u_1.dot( n_1 );
         un_2 = d2Length * u_2.dot( n_2 );

         if ( un_0 >= 0 )
         {
            c_up_0 = src[faceMemoryIndex];
         }
         else
         {
            // check if neighbor exists and get its id
            if ( adjEdges[0]->opposite_face_exists( face->getID() ) )
            {
               PrimitiveID oppositeFaceId          = adjEdges[0]->get_opposite_face( face->getID() );
               uint_t      localOppositeFaceId     = vertex.face_index( oppositeFaceId );
               uint_t      oppositeFaceMemoryIndex = 2 * localOppositeFaceId;
               c_up_0                              = src[oppositeFaceMemoryIndex];
            }
            else
            {
               // TODO: Handle boundary conditions in this case?
               c_up_0 = 0.0;
            }
         }

         if ( un_1 >= 0 )
         {
            c_up_1 = src[faceMemoryIndex];
         }
         else
         {
            // check if neighbor exists and get its id
            if ( adjEdges[1]->opposite_face_exists( face->getID() ) )
            {
               PrimitiveID oppositeFaceId          = adjEdges[1]->get_opposite_face( face->getID() );
               uint_t      localOppositeFaceId     = vertex.face_index( oppositeFaceId );
               uint_t      oppositeFaceMemoryIndex = 2 * localOppositeFaceId;
               c_up_1                              = src[oppositeFaceMemoryIndex];
            }
            else
            {
               // TODO: Handle boundary conditions in this case?
               c_up_1 = 0.0;
            }
         }

         if ( un_2 >= 0 )
         {
            c_up_2 = src[faceMemoryIndex];
         }
         else
         {
            c_up_2 = src[blueMemoryIndex];
         }

         tmp = un_0 * c_up_0 + un_1 * c_up_1 + un_2 * c_up_2;
         tmp *= faceAreaInv;

         if ( updateType == Replace )
         {
            dst[faceMemoryIndex] = tmp;
         }
         else if ( updateType == Add )
         {
            dst[faceMemoryIndex] += tmp;
         }
      }
   }
};

} // namespace hyteg
