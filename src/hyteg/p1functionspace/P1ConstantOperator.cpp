/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "P1ConstantOperator.hpp"

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#include "hyteg/fenics/fenics.hpp"
// #include "hyteg/forms/form_fenics_generated/p1_diffusion.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include "core/OpenMP.h"

#include "hyteg/forms/P1RowSumForm.hpp"
#include "hyteg/forms/form_fenics_base/P1ToP2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2ToP1FenicsForm.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_2D_macroface_vertexdof_to_vertexdof_add.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_2D_macroface_vertexdof_to_vertexdof_replace.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_3D_macrocell_vertexdof_to_vertexdof_add.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_3D_macrocell_vertexdof_to_vertexdof_replace.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace.hpp"
#include "hyteg/p1functionspace/generatedKernels/gaussseidel_3D_macrocell_P1.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_2D_macroface_vertexdof_to_vertexdof.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_2D_macroface_vertexdof_to_vertexdof_backwards.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_3D_macrocell_P1.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_3D_macrocell_P1_backwards.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_3D_macroface_P1.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_3D_macroface_P1_backwards.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_3D_macroface_P1_one_sided.hpp"
#include "hyteg/p1functionspace/generatedKernels/sor_3D_macroface_P1_one_sided_backwards.hpp"
#include "hyteg/p1functionspace/variablestencil/VertexDoFVariableStencil.hpp"

#include "P1Elements.hpp"

namespace hyteg {

using walberla::int_c;

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::P1ConstantOperator(
    const std::shared_ptr< PrimitiveStorage >& storage,
    size_t                                     minLevel,
    size_t                                     maxLevel )
: P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >( storage, minLevel, maxLevel, P1Form() )
{}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::P1ConstantOperator(
    const std::shared_ptr< PrimitiveStorage >& storage,
    size_t                                     minLevel,
    size_t                                     maxLevel,
    const P1Form&                              form )
: Operator( storage, minLevel, maxLevel )
, form_( form )
{
   auto cellP1StencilMemoryDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell > >(
           minLevel_, maxLevel_ );

   auto face3DP1StencilMemoryDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< vertexdof::macroface::StencilMap_T >, Face > >(
           minLevel_, maxLevel_ );

   auto edge3DP1StencilMemoryDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< vertexdof::macroedge::StencilMap_T >, Edge > >(
           minLevel_, maxLevel_ );

   auto faceP1StencilMemoryDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Face > >(
       minLevel_, maxLevel_, vertexDoFMacroFaceStencilMemorySize );
   auto edgeP1StencilMemoryDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Edge > >(
       minLevel_, maxLevel_, vertexDoFMacroEdgeStencilMemorySize );
   auto vertexP1StencilMemoryDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Vertex > >(
       minLevel_, maxLevel_, vertexDoFMacroVertexStencilMemorySize );

   storage->addCellData( cellStencilID_, cellP1StencilMemoryDataHandling, "P1OperatorCellStencil" );
   storage->addFaceData( faceStencilID_, faceP1StencilMemoryDataHandling, "P1OperatorFaceStencil" );
   storage->addFaceData( faceStencil3DID_, face3DP1StencilMemoryDataHandling, "P1OperatorFace3DStencil" );
   storage->addEdgeData( edgeStencilID_, edgeP1StencilMemoryDataHandling, "P1OperatorEdgeStencil" );
   storage->addEdgeData( edgeStencil3DID_, edge3DP1StencilMemoryDataHandling, "P1OperatorEdge3DStencil" );
   storage->addVertexData( vertexStencilID_, vertexP1StencilMemoryDataHandling, "P1OperatorVertexStencil" );

   if ( storage_->hasGlobalCells() )
   {
      const bool assemblyDefined = form_.assembly3DDefined();
      WALBERLA_CHECK( assemblyDefined, "Assembly undefined for 3D elements." );
      if ( form_.assemble3D() )
      {
         assembleStencils3D();
      }
   }
   else
   {
      if ( form_.assemble2D() )
      {
         const bool assemblyDefined = form_.assembly2DDefined();
         WALBERLA_CHECK( assemblyDefined, "Assembly undefined for 2D elements." );
         assembleStencils();
      }
   }
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::assembleStencils3D()
{
   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      for ( const auto& it : storage_->getVertices() )
      {
         auto vertex        = it.second;
         auto stencilSize   = vertex->getData( getVertexStencilID() )->getSize( level );
         auto stencilMemory = vertex->getData( getVertexStencilID() )->getPointer( level );

         form_.setGeometryMap( vertex->getGeometryMap() );
         auto stencil =
             P1Elements::P1Elements3D::assembleP1LocalStencil( storage_, *vertex, indexing::Index( 0, 0, 0 ), level, form_ );

         WALBERLA_ASSERT_EQUAL( stencilSize, stencil.size() );
         for ( uint_t i = 0; i < stencilSize; i++ )
         {
            stencilMemory[i] = stencil[i];
         }

         if ( Lumped )
         {
            for ( uint_t i = 1; i < stencilSize; i++ )
            {
               stencilMemory[0] += stencilMemory[i];
               stencilMemory[i] = 0;
            }
         }
         if ( Diagonal )
         {
            for ( uint_t i = 1; i < stencilSize; i++ )
            {
               stencilMemory[i] = 0;
            }
         }
         if ( InvertDiagonal )
         {
            stencilMemory[0] = 1.0 / stencilMemory[0];
         }
      }

      if ( level >= 1 )
      {
         for ( const auto& it : storage_->getEdges() )
         {
            auto  edge          = it.second;
            auto  stencilSize   = edge->getData( getEdgeStencilID() )->getSize( level );
            auto  stencilMemory = edge->getData( getEdgeStencilID() )->getPointer( level );
            auto& stencilMap    = edge->getData( getEdgeStencil3DID() )->getData( level );

            form_.setGeometryMap( edge->getGeometryMap() );

            // old linear stencil memory
            auto stencil =
                P1Elements::P1Elements3D::assembleP1LocalStencil( storage_, *edge, indexing::Index( 1, 0, 0 ), level, form_ );

            WALBERLA_ASSERT_EQUAL( stencilSize, stencil.size() );
            for ( uint_t i = 0; i < stencilSize; i++ )
            {
               stencilMemory[i] = stencil[i];
            }

            if ( Lumped )
            {
               stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )] +=
                   stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_W )];
               stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_W )] = 0;
               stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )] +=
                   stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_E )];
               stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_E )] = 0;
               for ( uint_t neighborFace = 0; neighborFace < it.second->getNumNeighborFaces(); neighborFace++ )
               {
                  stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )] +=
                      stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborFace( stencilDirection::VERTEX_W, neighborFace )];
                  stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborFace( stencilDirection::VERTEX_W, neighborFace )] = 0;
                  stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )] +=
                      stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborFace( stencilDirection::VERTEX_E, neighborFace )];
                  stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborFace( stencilDirection::VERTEX_E, neighborFace )] = 0;
               }
               for ( uint_t neighborCell = 0; neighborCell < it.second->getNumNeighborCells(); neighborCell++ )
               {
                  stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )] +=
                      stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell,
                                                                                      it.second->getNumNeighborFaces() )];
                  stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell,
                                                                                  it.second->getNumNeighborFaces() )] = 0;
               }
            }
            if ( Diagonal )
            {
               stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_W )] = 0;
               stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_E )] = 0;
               for ( uint_t neighborFace = 0; neighborFace < it.second->getNumNeighborFaces(); neighborFace++ )
               {
                  stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborFace( stencilDirection::VERTEX_W, neighborFace )] = 0;
                  stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborFace( stencilDirection::VERTEX_E, neighborFace )] = 0;
               }
               for ( uint_t neighborCell = 0; neighborCell < it.second->getNumNeighborCells(); neighborCell++ )
               {
                  stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell,
                                                                                  it.second->getNumNeighborFaces() )] = 0;
               }
            }
            if ( InvertDiagonal )
            {
               stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )] =
                   1.0 / stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )];
            }

            // new map
            for ( uint_t neighborCellID = 0; neighborCellID < edge->getNumNeighborCells(); neighborCellID++ )
            {
               auto neighborCell = storage_->getCell( edge->neighborCells().at( neighborCellID ) );
               auto vertexAssemblyIndexInCell =
                   vertexdof::macroedge::getIndexInNeighboringMacroCell( { 1, 0, 0 }, *edge, neighborCellID, *storage_, level );
               stencilMap[neighborCellID] = P1Elements::P1Elements3D::assembleP1LocalStencilNew(
                   storage_, *neighborCell, vertexAssemblyIndexInCell, level, form_ );
            }
         }
      }

      if ( level >= 1 )
      {
         for ( const auto& it : storage_->getFaces() )
         {
            auto  face          = it.second;
            auto& stencilMemory = face->getData( getFaceStencil3DID() )->getData( level );

            for ( uint_t neighborCellID = 0; neighborCellID < face->getNumNeighborCells(); neighborCellID++ )
            {
               auto neighborCell = storage_->getCell( face->neighborCells().at( neighborCellID ) );
               auto vertexAssemblyIndexInCell =
                   vertexdof::macroface::getIndexInNeighboringMacroCell( { 1, 1, 0 }, *face, neighborCellID, *storage_, level );
               stencilMemory[neighborCellID] = P1Elements::P1Elements3D::assembleP1LocalStencilNew(
                   storage_, *neighborCell, vertexAssemblyIndexInCell, level, form_ );
            }

            // The lumping and inverted diagonal modifications for split stencils is realized
            // by adding up the diagonal entries on the _first_ neighbor cell and setting all other parts to zero.

            WALBERLA_ASSERT_GREATER( face->getNumNeighborCells(), 0 );

            if ( Lumped )
            {
               for ( uint_t neighborCellID = 0; neighborCellID < face->getNumNeighborCells(); neighborCellID++ )
               {
                  for ( auto& stencilIt : stencilMemory[neighborCellID] )
                  {
                     if ( !( neighborCellID == 0 && stencilIt.first == indexing::IndexIncrement( { 0, 0, 0 } ) ) )
                     {
                        stencilMemory[0][{ 0, 0, 0 }] += stencilIt.second;
                        stencilIt.second = 0;
                     }
                  }
               }
            }
            if ( Diagonal )
            {
               for ( uint_t neighborCellID = 0; neighborCellID < face->getNumNeighborCells(); neighborCellID++ )
               {
                  for ( auto& stencilIt : stencilMemory[neighborCellID] )
                  {
                     if ( stencilIt.first != indexing::IndexIncrement( { 0, 0, 0 } ) )
                     {
                        stencilIt.second = 0;
                     }
                  }
               }
            }
            if ( InvertDiagonal )
            {
               for ( uint_t neighborCellID = 1; neighborCellID < face->getNumNeighborCells(); neighborCellID++ )
               {
                  stencilMemory[0][{ 0, 0, 0 }] += stencilMemory[neighborCellID][{ 0, 0, 0 }];
                  stencilMemory[neighborCellID][{ 0, 0, 0 }] = 0;
               }
               stencilMemory[0][{ 0, 0, 0 }] = 1.0 / stencilMemory[0][{ 0, 0, 0 }];
            }
         }
      }

      if ( level >= 2 )
      {
         for ( const auto& it : storage_->getCells() )
         {
            auto  cell          = it.second;
            auto& stencilMemory = cell->getData( getCellStencilID() )->getData( level );

            form_.setGeometryMap( cell->getGeometryMap() );

            stencilMemory =
                P1Elements::P1Elements3D::assembleP1LocalStencilNew( storage_, *cell, indexing::Index( 1, 1, 1 ), level, form_ );

            if ( Lumped )
            {
               for ( auto dir : vertexdof::macrocell::neighborsWithoutCenter )
               {
                  stencilMemory[{ 0, 0, 0 }] += stencilMemory[vertexdof::logicalIndexOffsetFromVertex( dir )];
                  stencilMemory[vertexdof::logicalIndexOffsetFromVertex( dir )] = 0;
               }
            }
            if ( Diagonal )
            {
               for ( auto dir : vertexdof::macrocell::neighborsWithoutCenter )
               {
                  stencilMemory[vertexdof::logicalIndexOffsetFromVertex( dir )] = 0;
               }
            }
            if ( InvertDiagonal )
            {
               stencilMemory[{ 0, 0, 0 }] = 1.0 / stencilMemory[{ 0, 0, 0 }];
            }
         }
      }
   }
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::assembleStencils()
{
   using namespace P1Elements::P1Elements2D;
   typedef stencilDirection sD;

   for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
   {
      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         auto face_stencil = face.getData( faceStencilID_ )->getPointer( level );

         Point3D x( face.coords[0] );
         uint_t  rowsize = levelinfo::num_microvertices_per_edge( level );
         real_t  h       = 1.0 / ( walberla::real_c( rowsize - 1 ) );

         Point3D d0 = h * ( face.coords[1] - face.coords[0] );
         Point3D d2 = h * ( face.coords[2] - face.coords[0] );

         form_.setGeometryMap( face.getGeometryMap() );

         Point3D dirS  = -1.0 * d2;
         Point3D dirSE = d0 - 1.0 * d2;
         Point3D dirE  = d0;
         Point3D dirW  = -1.0 * d0;
         Point3D dirNW = -1.0 * d0 + d2;
         Point3D dirN  = d2;

         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form_, { x, x + dirW, x + dirS }, P1Elements::P1Elements2D::elementSW, face_stencil );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form_, { x, x + dirS, x + dirSE }, P1Elements::P1Elements2D::elementS, face_stencil );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form_, { x, x + dirSE, x + dirE }, P1Elements::P1Elements2D::elementSE, face_stencil );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form_, { x, x + dirE, x + dirN }, P1Elements::P1Elements2D::elementNE, face_stencil );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form_, { x, x + dirN, x + dirNW }, P1Elements::P1Elements2D::elementN, face_stencil );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form_, { x, x + dirNW, x + dirW }, P1Elements::P1Elements2D::elementNW, face_stencil );

         if ( Lumped )
         {
            for ( const auto& neighbor : vertexdof::macroface::neighborsWithoutCenter )
            {
               face_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] +=
                   face_stencil[vertexdof::stencilIndexFromVertex( neighbor )];
               face_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
            }
         }

         if ( Diagonal )
         {
            for ( const auto& neighbor : vertexdof::macroface::neighborsWithoutCenter )
            {
               face_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
            }
         }

         if ( InvertDiagonal )
         {
            face_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] =
                1.0 / face_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )];
         }
      }

      for ( auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         auto edge_stencil = edge.getData( edgeStencilID_ )->getPointer( level );

         size_t rowsize = levelinfo::num_microvertices_per_edge( level );

         Face*  faceS   = storage_->getFace( edge.neighborFaces()[0] );
         Face*  faceN   = nullptr;
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

         Point3D x = edge.getCoordinates()[0];

         uint_t  s_north, e_north, o_north;
         Point3D dir_N;
         Point3D dir_NW;

         if ( edge.getNumNeighborFaces() == 2 )
         {
            faceN   = storage_->getFace( edge.neighborFaces()[1] );
            s_north = faceN->vertex_index( edge.neighborVertices()[0] );
            e_north = faceN->vertex_index( edge.neighborVertices()[1] );
            o_north = faceN->vertex_index( faceN->get_vertex_opposite_to_edge( edge.getID() ) );

            Point3D dN_so = h * ( faceN->coords[o_north] - faceN->coords[s_north] );
            Point3D dN_oe = h * ( faceN->coords[e_north] - faceN->coords[o_north] );

            dir_N  = dN_so;
            dir_NW = -1.0 * dN_oe;
         }

         // assemble south
         form_.setGeometryMap( faceS->getGeometryMap() );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form_, { x, x + dir_W, x + dir_S }, P1Elements::P1Elements2D::elementSW, edge_stencil );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form_, { x, x + dir_S, x + dir_SE }, P1Elements::P1Elements2D::elementS, edge_stencil );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form_, { x, x + dir_SE, x + dir_E }, P1Elements::P1Elements2D::elementSE, edge_stencil );

         if ( edge.getNumNeighborFaces() == 2 )
         {
            form_.setGeometryMap( faceN->getGeometryMap() );
            vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                form_, { x, x + dir_E, x + dir_N }, P1Elements::P1Elements2D::elementNE, edge_stencil );
            vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                form_, { x, x + dir_N, x + dir_NW }, P1Elements::P1Elements2D::elementN, edge_stencil );
            vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                form_, { x, x + dir_NW, x + dir_W }, P1Elements::P1Elements2D::elementNW, edge_stencil );
         }

         if ( Lumped )
         {
            for ( const auto& neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
            {
               edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] +=
                   edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )];
               edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
            }

            for ( const auto& neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
            {
               edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] +=
                   edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )];
               edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
            }

            if ( edge.getNumNeighborFaces() == 2 )
            {
               for ( const auto& neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
               {
                  edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] +=
                      edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )];
                  edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
               }
            }
         }

         if ( Diagonal )
         {
            for ( const auto& neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
            {
               edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
            }

            for ( const auto& neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
            {
               edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
            }

            if ( edge.getNumNeighborFaces() == 2 )
            {
               for ( const auto& neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
               {
                  edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
               }
            }
         }

         if ( InvertDiagonal )
         {
            edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] =
                1.0 / edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )];
         }
      }

      for ( auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         auto vertex_stencil = vertex.getData( vertexStencilID_ )->getPointer( level );

         uint_t rowsize = levelinfo::num_microvertices_per_edge( level );

         Point3D x;
         Point3D d0;
         Point3D d2;

         real_t h = 1.0 / ( walberla::real_c( rowsize - 1 ) );

         uint_t neighborId = 0;
         for ( auto& faceId : vertex.neighborFaces() )
         {
            Face* face = storage_->getFace( faceId );
            form_.setGeometryMap( face->getGeometryMap() );

            uint_t                     v_i       = face->vertex_index( vertex.getID() );
            std::vector< PrimitiveID > adj_edges = face->adjacent_edges( vertex.getID() );

            x  = face->coords[v_i];
            d0 = ( face->coords[face->vertex_index( storage_->getEdge( adj_edges[0] )->get_opposite_vertex( vertex.getID() ) )] -
                   x ) *
                 h;
            d2 = ( face->coords[face->vertex_index( storage_->getEdge( adj_edges[1] )->get_opposite_vertex( vertex.getID() ) )] -
                   x ) *
                 h;

            Point3D matrixRow;
            form_.integrate( { { x, x + d0, x + d2 } }, matrixRow );

            uint_t i = 1;
            // iterate over adjacent edges
            for ( auto& edgeId : adj_edges )
            {
               uint_t edge_idx = vertex.edge_index( edgeId ) + 1;

               vertex_stencil[edge_idx] += matrixRow[i];
               i += 1;
            }

            // add contribution of center vertex
            vertex_stencil[0] += matrixRow[0];

            ++neighborId;
         }

         if ( Lumped )
         {
            for ( uint_t i = 1; i < vertex.getData( vertexStencilID_ )->getSize( level ); ++i )
            {
               vertex_stencil[0] += vertex_stencil[i];
               vertex_stencil[i] = 0;
            }
         }

         if ( Diagonal )
         {
            for ( uint_t i = 1; i < vertex.getData( vertexStencilID_ )->getSize( level ); ++i )
            {
               vertex_stencil[i] = 0;
            }
         }

         if ( InvertDiagonal )
         {
            vertex_stencil[0] = 1.0 / vertex_stencil[0];
         }
      }
   }
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::apply( const P1Function< real_t >& src,
                                                                            const P1Function< real_t >& dst,
                                                                            size_t                      level,
                                                                            DoFType                     flag,
                                                                            UpdateType                  updateType ) const
{
   WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

   this->startTiming( "Apply" );
   src.communicate< Vertex, Edge >( level );
   src.communicate< Edge, Face >( level );
   src.communicate< Face, Cell >( level );

   src.communicate< Cell, Face >( level );
   src.communicate< Face, Edge >( level );
   src.communicate< Edge, Vertex >( level );

   this->timingTree_->start( "Macro-Vertex" );

   std::vector< PrimitiveID > vertexIDs = this->getStorage()->getVertexIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( vertexIDs.size() ); i++ )
   {
      Vertex& vertex = *this->getStorage()->getVertex( vertexIDs[uint_c( i )] );

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         vertexdof::macrovertex::apply< real_t >(
             vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), level, updateType );
      }
   }

   this->timingTree_->stop( "Macro-Vertex" );

   this->timingTree_->start( "Macro-Edge" );

   if ( level >= 1 )
   {
      std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
      for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
      {
         Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            vertexdof::macroedge::apply< real_t >(
                level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
         }
      }
   }

   this->timingTree_->stop( "Macro-Edge" );

   this->timingTree_->start( "Macro-Face" );

   if ( level >= 2 )
   {
      std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
      for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
      {
         Face& face = *this->getStorage()->getFace( faceIDs[uint_c( i )] );

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            if ( storage_->hasGlobalCells() )
            {
               if ( hyteg::globalDefines::useGeneratedKernels )
               {
                  if ( face.getNumNeighborCells() == 2 )
                  {
                     WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start( "Two-sided" ); }
                  }
                  else
                  {
                     WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start( "One-sided" ); }
                  }

                  auto         opr_data    = face.getData( faceStencil3DID_ )->getData( level );
                  auto         src_data    = face.getData( src.getFaceDataID() )->getPointer( level );
                  auto         dst_data    = face.getData( dst.getFaceDataID() )->getPointer( level );
                  const uint_t offset_gl_0 = levelinfo::num_microvertices_per_face( level );

                  auto neighborCell0 = storage_->getCell( face.neighborCells()[0] );

                  auto neighbor_cell_0_local_vertex_id_0 =
                      static_cast< int32_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                                  .at( neighborCell0->getLocalFaceID( face.getID() ) )
                                                  .at( 0 ) );
                  auto neighbor_cell_0_local_vertex_id_1 =
                      static_cast< int32_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                                  .at( neighborCell0->getLocalFaceID( face.getID() ) )
                                                  .at( 1 ) );
                  auto neighbor_cell_0_local_vertex_id_2 =
                      static_cast< int32_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                                  .at( neighborCell0->getLocalFaceID( face.getID() ) )
                                                  .at( 2 ) );

                  if ( updateType == Replace )
                  {
                     vertexdof::macroface::generated::apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace(
                         dst_data,
                         src_data,
                         &src_data[offset_gl_0],
                         static_cast< int32_t >( level ),
                         neighbor_cell_0_local_vertex_id_0,
                         neighbor_cell_0_local_vertex_id_1,
                         neighbor_cell_0_local_vertex_id_2,
                         opr_data[0] );
                  }
                  else
                  {
                     vertexdof::macroface::generated::apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add(
                         dst_data,
                         src_data,
                         &src_data[offset_gl_0],
                         static_cast< int32_t >( level ),
                         neighbor_cell_0_local_vertex_id_0,
                         neighbor_cell_0_local_vertex_id_1,
                         neighbor_cell_0_local_vertex_id_2,
                         opr_data[0] );
                  }

                  if ( face.getNumNeighborCells() == 2 )
                  {
                     const uint_t offset_gl_1 = offset_gl_0 + levelinfo::num_microvertices_per_face_from_width(
                                                                  levelinfo::num_microvertices_per_edge( level ) - 1 );

                     auto neighborCell1 = storage_->getCell( face.neighborCells()[1] );

                     auto neighbor_cell_1_local_vertex_id_0 =
                         static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                                     .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                                     .at( 0 ) );
                     auto neighbor_cell_1_local_vertex_id_1 =
                         static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                                     .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                                     .at( 1 ) );
                     auto neighbor_cell_1_local_vertex_id_2 =
                         static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                                     .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                                     .at( 2 ) );

                     vertexdof::macroface::generated::apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add(
                         dst_data,
                         src_data,
                         &src_data[offset_gl_1],
                         static_cast< int32_t >( level ),
                         neighbor_cell_1_local_vertex_id_0,
                         neighbor_cell_1_local_vertex_id_1,
                         neighbor_cell_1_local_vertex_id_2,
                         opr_data[1] );
                  }

                  if ( face.getNumNeighborCells() == 2 )
                  {
                     WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->stop( "Two-sided" ); }
                  }
                  else
                  {
                     WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->stop( "One-sided" ); }
                  }
               }
               else
               {
                  vertexdof::macroface::apply3D< real_t >(
                      level, face, *storage_, faceStencil3DID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
               }
            }
            else
            {
               if ( hyteg::globalDefines::useGeneratedKernels )
               {
                  real_t* opr_data = face.getData( faceStencilID_ )->getPointer( level );
                  real_t* src_data = face.getData( src.getFaceDataID() )->getPointer( level );
                  real_t* dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );
                  if ( updateType == hyteg::Replace )
                  {
                     vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_replace(
                         dst_data, src_data, opr_data, static_cast< int32_t >( level ) );
                  }
                  else if ( updateType == hyteg::Add )
                  {
                     vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_add(
                         dst_data, src_data, opr_data, static_cast< int32_t >( level ) );
                  }
               }
               else
               {
                  vertexdof::macroface::apply< real_t >(
                      level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
               }
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Face" );

   this->timingTree_->start( "Macro-Cell" );

   if ( level >= 2 )
   {
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
      for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
      {
         Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );

         const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
         if ( testFlag( cellBC, flag ) )
         {
            if ( hyteg::globalDefines::useGeneratedKernels )
            {
               auto    opr_data = cell.getData( cellStencilID_ )->getData( level );
               real_t* src_data = cell.getData( src.getCellDataID() )->getPointer( level );
               real_t* dst_data = cell.getData( dst.getCellDataID() )->getPointer( level );
               if ( updateType == Replace )
               {
                  vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_replace(
                      dst_data, src_data, static_cast< int32_t >( level ), opr_data );
               }
               else if ( updateType == Add )
               {
                  vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_add(
                      dst_data, src_data, static_cast< int32_t >( level ), opr_data );
               }
            }
            else
            {
               vertexdof::macrocell::apply< real_t >(
                   level, cell, cellStencilID_, src.getCellDataID(), dst.getCellDataID(), updateType );
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Cell" );

   this->stopTiming( "Apply" );
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::smooth_gs( const P1Function< real_t >& dst,
                                                                                const P1Function< real_t >& rhs,
                                                                                size_t                      level,
                                                                                DoFType                     flag ) const
{
   this->startTiming( "Gauss-Seidel" );

   dst.communicate< Vertex, Edge >( level );
   dst.communicate< Edge, Face >( level );
   dst.communicate< Face, Cell >( level );

   dst.communicate< Cell, Face >( level );
   dst.communicate< Face, Edge >( level );
   dst.communicate< Edge, Vertex >( level );

   this->timingTree_->start( "Macro-Vertex" );

   for ( auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         vertexdof::macrovertex::smooth_sor( vertex, vertexStencilID_, dst.getVertexDataID(), rhs.getVertexDataID(), level, 1.0 );
      }
   }

   this->timingTree_->stop( "Macro-Vertex" );

   dst.communicate< Vertex, Edge >( level );

   this->timingTree_->start( "Macro-Edge" );

   for ( auto& it : storage_->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         vertexdof::macroedge::smooth_sor< real_t >( level, edge, edgeStencilID_, dst.getEdgeDataID(), rhs.getEdgeDataID(), 1.0 );
      }
   }

   this->timingTree_->stop( "Macro-Edge" );

   dst.communicate< Edge, Face >( level );

   this->timingTree_->start( "Macro-Face" );

   for ( auto& it : storage_->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         if ( storage_->hasGlobalCells() )
         {
            if ( globalDefines::useGeneratedKernels && face.getNumNeighborCells() == 2 )
            {
               this->timingTree_->start( "Two-sided" );
               auto rhs_data = face.getData( rhs.getFaceDataID() )->getPointer( level );
               auto dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );
               auto stencil  = face.getData( faceStencil3DID_ )->getData( level );

               auto neighborCell0 = storage_->getCell( face.neighborCells()[0] );
               auto neighborCell1 = storage_->getCell( face.neighborCells()[1] );

               auto neighbor_cell_0_local_vertex_id_0 =
                   static_cast< int32_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                               .at( neighborCell0->getLocalFaceID( face.getID() ) )
                                               .at( 0 ) );
               auto neighbor_cell_0_local_vertex_id_1 =
                   static_cast< int32_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                               .at( neighborCell0->getLocalFaceID( face.getID() ) )
                                               .at( 1 ) );
               auto neighbor_cell_0_local_vertex_id_2 =
                   static_cast< int32_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                               .at( neighborCell0->getLocalFaceID( face.getID() ) )
                                               .at( 2 ) );

               auto neighbor_cell_1_local_vertex_id_0 =
                   static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                               .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                               .at( 0 ) );
               auto neighbor_cell_1_local_vertex_id_1 =
                   static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                               .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                               .at( 1 ) );
               auto neighbor_cell_1_local_vertex_id_2 =
                   static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                               .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                               .at( 2 ) );

               const uint_t vertex_offset_gl_0 = levelinfo::num_microvertices_per_face( level );
               const uint_t vertex_offset_gl_1 = vertex_offset_gl_0 + levelinfo::num_microvertices_per_face_from_width(
                                                                          levelinfo::num_microvertices_per_edge( level ) - 1 );

               if ( neighbor_cell_0_local_vertex_id_0 > neighbor_cell_1_local_vertex_id_0 ||
                    ( neighbor_cell_0_local_vertex_id_0 == neighbor_cell_1_local_vertex_id_0 &&
                      neighbor_cell_0_local_vertex_id_1 > neighbor_cell_1_local_vertex_id_1 ) ||
                    ( neighbor_cell_0_local_vertex_id_0 == neighbor_cell_1_local_vertex_id_0 &&
                      neighbor_cell_0_local_vertex_id_1 == neighbor_cell_1_local_vertex_id_1 &&
                      neighbor_cell_0_local_vertex_id_2 > neighbor_cell_1_local_vertex_id_2 ) )
               {
                  vertexdof::macroface::generated::sor_3D_macroface_P1( dst_data,
                                                                        &dst_data[vertex_offset_gl_1],
                                                                        &dst_data[vertex_offset_gl_0],
                                                                        rhs_data,
                                                                        static_cast< int32_t >( level ),
                                                                        neighbor_cell_1_local_vertex_id_0,
                                                                        neighbor_cell_1_local_vertex_id_1,
                                                                        neighbor_cell_1_local_vertex_id_2,
                                                                        neighbor_cell_0_local_vertex_id_0,
                                                                        neighbor_cell_0_local_vertex_id_1,
                                                                        neighbor_cell_0_local_vertex_id_2,
                                                                        1.0,
                                                                        stencil[1],
                                                                        stencil[0] );
               }
               else
               {
                  vertexdof::macroface::generated::sor_3D_macroface_P1( dst_data,
                                                                        &dst_data[vertex_offset_gl_0],
                                                                        &dst_data[vertex_offset_gl_1],
                                                                        rhs_data,
                                                                        static_cast< int32_t >( level ),
                                                                        neighbor_cell_0_local_vertex_id_0,
                                                                        neighbor_cell_0_local_vertex_id_1,
                                                                        neighbor_cell_0_local_vertex_id_2,
                                                                        neighbor_cell_1_local_vertex_id_0,
                                                                        neighbor_cell_1_local_vertex_id_1,
                                                                        neighbor_cell_1_local_vertex_id_2,
                                                                        1.0,
                                                                        stencil[0],
                                                                        stencil[1] );
               }
               this->timingTree_->stop( "Two-sided" );
            }
            else
            {
               this->timingTree_->start( "One-sided" );
               vertexdof::macroface::smoothSOR3D< real_t >(
                   level, face, *storage_, faceStencil3DID_, dst.getFaceDataID(), rhs.getFaceDataID(), 1.0 );
               this->timingTree_->stop( "One-sided" );
            }
         }
         else
         {
            vertexdof::macroface::smooth_sor< real_t >(
                level, face, faceStencilID_, dst.getFaceDataID(), rhs.getFaceDataID(), 1.0 );
         }
      }
   }

   this->timingTree_->stop( "Macro-Face" );

   dst.communicate< Face, Cell >( level );

   this->timingTree_->start( "Macro-Cell" );

   for ( auto& it : storage_->getCells() )
   {
      Cell& cell = *it.second;

      const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
      if ( testFlag( cellBC, flag ) )
      {
         if ( globalDefines::useGeneratedKernels )
         {
            auto rhs_data = cell.getData( rhs.getCellDataID() )->getPointer( level );
            auto dst_data = cell.getData( dst.getCellDataID() )->getPointer( level );
            auto stencil  = cell.getData( cellStencilID_ )->getData( level );
            vertexdof::macrocell::generated::gaussseidel_3D_macrocell_P1(
                dst_data, rhs_data, static_cast< int32_t >( level ), stencil );
         }
         else
         {
            vertexdof::macrocell::smooth_sor< real_t >(
                level, cell, cellStencilID_, dst.getCellDataID(), rhs.getCellDataID(), 1.0 );
         }
      }
   }

   this->timingTree_->stop( "Macro-Cell" );

   this->stopTiming( "Gauss-Seidel" );
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::smooth_sor_macro_vertices( const P1Function< real_t >& dst,
                                                                                                const P1Function< real_t >& rhs,
                                                                                                real_t                      relax,
                                                                                                size_t                      level,
                                                                                                DoFType                     flag,
                                                                                                const bool& backwards ) const
{
   WALBERLA_UNUSED( backwards );

   this->timingTree_->start( "Macro-Vertex" );

   for ( auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         vertexdof::macrovertex::smooth_sor(
             vertex, vertexStencilID_, dst.getVertexDataID(), rhs.getVertexDataID(), level, relax );
      }
   }

   this->timingTree_->stop( "Macro-Vertex" );
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::smooth_sor_macro_edges( const P1Function< real_t >& dst,
                                                                                             const P1Function< real_t >& rhs,
                                                                                             real_t                      relax,
                                                                                             size_t                      level,
                                                                                             DoFType                     flag,
                                                                                             const bool& backwards ) const
{
   this->timingTree_->start( "Macro-Edge" );

   for ( auto& it : storage_->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         vertexdof::macroedge::smooth_sor< real_t >(
             level, edge, edgeStencilID_, dst.getEdgeDataID(), rhs.getEdgeDataID(), relax, backwards );
      }
   }

   this->timingTree_->stop( "Macro-Edge" );
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::smooth_sor_macro_faces( const P1Function< real_t >& dst,
                                                                                             const P1Function< real_t >& rhs,
                                                                                             real_t                      relax,
                                                                                             size_t                      level,
                                                                                             DoFType                     flag,
                                                                                             const bool& backwards ) const
{
   this->timingTree_->start( "Macro-Face" );

   for ( auto& it : storage_->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         if ( storage_->hasGlobalCells() )
         {
            if ( globalDefines::useGeneratedKernels )
            {
               auto rhs_data = face.getData( rhs.getFaceDataID() )->getPointer( level );
               auto dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );
               auto stencil  = face.getData( faceStencil3DID_ )->getData( level );

               auto neighborCell0 = storage_->getCell( face.neighborCells()[0] );

               auto neighbor_cell_0_local_vertex_id_0 =
                   static_cast< int32_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                               .at( neighborCell0->getLocalFaceID( face.getID() ) )
                                               .at( 0 ) );
               auto neighbor_cell_0_local_vertex_id_1 =
                   static_cast< int32_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                               .at( neighborCell0->getLocalFaceID( face.getID() ) )
                                               .at( 1 ) );
               auto neighbor_cell_0_local_vertex_id_2 =
                   static_cast< int32_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps()
                                               .at( neighborCell0->getLocalFaceID( face.getID() ) )
                                               .at( 2 ) );

               const uint_t vertex_offset_gl_0 = levelinfo::num_microvertices_per_face( level );

               if ( face.getNumNeighborCells() == 1 )
               {
                  this->timingTree_->start( "One-sided" );
                  if ( backwards )
                  {
                     vertexdof::macroface::generated::sor_3D_macroface_P1_one_sided_backwards( dst_data,
                                                                                               &dst_data[vertex_offset_gl_0],
                                                                                               rhs_data,
                                                                                               static_cast< int32_t >( level ),
                                                                                               neighbor_cell_0_local_vertex_id_0,
                                                                                               neighbor_cell_0_local_vertex_id_1,
                                                                                               neighbor_cell_0_local_vertex_id_2,
                                                                                               relax,
                                                                                               stencil[0] );
                  }
                  else
                  {
                     vertexdof::macroface::generated::sor_3D_macroface_P1_one_sided( dst_data,
                                                                                     &dst_data[vertex_offset_gl_0],
                                                                                     rhs_data,
                                                                                     static_cast< int32_t >( level ),
                                                                                     neighbor_cell_0_local_vertex_id_0,
                                                                                     neighbor_cell_0_local_vertex_id_1,
                                                                                     neighbor_cell_0_local_vertex_id_2,
                                                                                     relax,
                                                                                     stencil[0] );
                  }
                  this->timingTree_->stop( "One-sided" );
               }
               if ( face.getNumNeighborCells() == 2 )
               {
                  this->timingTree_->start( "Two-sided" );

                  auto neighborCell1 = storage_->getCell( face.neighborCells()[1] );

                  auto neighbor_cell_1_local_vertex_id_0 =
                      static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                                  .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                                  .at( 0 ) );
                  auto neighbor_cell_1_local_vertex_id_1 =
                      static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                                  .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                                  .at( 1 ) );
                  auto neighbor_cell_1_local_vertex_id_2 =
                      static_cast< int32_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps()
                                                  .at( neighborCell1->getLocalFaceID( face.getID() ) )
                                                  .at( 2 ) );
                  const uint_t vertex_offset_gl_1 = vertex_offset_gl_0 + levelinfo::num_microvertices_per_face_from_width(
                                                                             levelinfo::num_microvertices_per_edge( level ) - 1 );

                  if ( neighbor_cell_0_local_vertex_id_0 > neighbor_cell_1_local_vertex_id_0 ||
                       ( neighbor_cell_0_local_vertex_id_0 == neighbor_cell_1_local_vertex_id_0 &&
                         neighbor_cell_0_local_vertex_id_1 > neighbor_cell_1_local_vertex_id_1 ) ||
                       ( neighbor_cell_0_local_vertex_id_0 == neighbor_cell_1_local_vertex_id_0 &&
                         neighbor_cell_0_local_vertex_id_1 == neighbor_cell_1_local_vertex_id_1 &&
                         neighbor_cell_0_local_vertex_id_2 > neighbor_cell_1_local_vertex_id_2 ) )
                  {
                     if ( backwards )
                     {
                        vertexdof::macroface::generated::sor_3D_macroface_P1_backwards( dst_data,
                                                                                        &dst_data[vertex_offset_gl_1],
                                                                                        &dst_data[vertex_offset_gl_0],
                                                                                        rhs_data,
                                                                                        static_cast< int32_t >( level ),
                                                                                        neighbor_cell_1_local_vertex_id_0,
                                                                                        neighbor_cell_1_local_vertex_id_1,
                                                                                        neighbor_cell_1_local_vertex_id_2,
                                                                                        neighbor_cell_0_local_vertex_id_0,
                                                                                        neighbor_cell_0_local_vertex_id_1,
                                                                                        neighbor_cell_0_local_vertex_id_2,
                                                                                        relax,
                                                                                        stencil[1],
                                                                                        stencil[0] );
                     }
                     else
                     {
                        vertexdof::macroface::generated::sor_3D_macroface_P1( dst_data,
                                                                              &dst_data[vertex_offset_gl_1],
                                                                              &dst_data[vertex_offset_gl_0],
                                                                              rhs_data,
                                                                              static_cast< int32_t >( level ),
                                                                              neighbor_cell_1_local_vertex_id_0,
                                                                              neighbor_cell_1_local_vertex_id_1,
                                                                              neighbor_cell_1_local_vertex_id_2,
                                                                              neighbor_cell_0_local_vertex_id_0,
                                                                              neighbor_cell_0_local_vertex_id_1,
                                                                              neighbor_cell_0_local_vertex_id_2,
                                                                              relax,
                                                                              stencil[1],
                                                                              stencil[0] );
                     }
                  }
                  else
                  {
                     if ( backwards )
                     {
                        vertexdof::macroface::generated::sor_3D_macroface_P1_backwards( dst_data,
                                                                                        &dst_data[vertex_offset_gl_0],
                                                                                        &dst_data[vertex_offset_gl_1],
                                                                                        rhs_data,
                                                                                        static_cast< int32_t >( level ),
                                                                                        neighbor_cell_0_local_vertex_id_0,
                                                                                        neighbor_cell_0_local_vertex_id_1,
                                                                                        neighbor_cell_0_local_vertex_id_2,
                                                                                        neighbor_cell_1_local_vertex_id_0,
                                                                                        neighbor_cell_1_local_vertex_id_1,
                                                                                        neighbor_cell_1_local_vertex_id_2,
                                                                                        relax,
                                                                                        stencil[0],
                                                                                        stencil[1] );
                     }
                     else
                     {
                        vertexdof::macroface::generated::sor_3D_macroface_P1( dst_data,
                                                                              &dst_data[vertex_offset_gl_0],
                                                                              &dst_data[vertex_offset_gl_1],
                                                                              rhs_data,
                                                                              static_cast< int32_t >( level ),
                                                                              neighbor_cell_0_local_vertex_id_0,
                                                                              neighbor_cell_0_local_vertex_id_1,
                                                                              neighbor_cell_0_local_vertex_id_2,
                                                                              neighbor_cell_1_local_vertex_id_0,
                                                                              neighbor_cell_1_local_vertex_id_1,
                                                                              neighbor_cell_1_local_vertex_id_2,
                                                                              relax,
                                                                              stencil[0],
                                                                              stencil[1] );
                     }
                  }
                  this->timingTree_->stop( "Two-sided" );
               }
            }
            else
            {
               vertexdof::macroface::smoothSOR3D< real_t >(
                   level, face, *storage_, faceStencil3DID_, dst.getFaceDataID(), rhs.getFaceDataID(), relax );
            }
         }
         else
         {
            if ( globalDefines::useGeneratedKernels )
            {
               auto rhs_data = face.getData( rhs.getFaceDataID() )->getPointer( level );
               auto dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );
               auto stencil  = face.getData( faceStencilID_ )->getPointer( level );

               if ( backwards )
               {
                  vertexdof::macroface::generated::sor_2D_macroface_vertexdof_to_vertexdof_backwards(
                      dst_data, rhs_data, stencil, static_cast< int32_t >( level ), relax );
               }
               else
               {
                  vertexdof::macroface::generated::sor_2D_macroface_vertexdof_to_vertexdof(
                      dst_data, rhs_data, stencil, static_cast< int32_t >( level ), relax );
               }
            }
            else
            {
               vertexdof::macroface::smooth_sor< real_t >(
                   level, face, faceStencilID_, dst.getFaceDataID(), rhs.getFaceDataID(), relax );
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Face" );
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::smooth_sor_macro_cells( const P1Function< real_t >& dst,
                                                                                             const P1Function< real_t >& rhs,
                                                                                             real_t                      relax,
                                                                                             size_t                      level,
                                                                                             DoFType                     flag,
                                                                                             const bool& backwards ) const
{
   this->timingTree_->start( "Macro-Cell" );

   for ( auto& it : storage_->getCells() )
   {
      Cell& cell = *it.second;

      const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
      if ( testFlag( cellBC, flag ) )
      {
         if ( globalDefines::useGeneratedKernels )
         {
            auto rhs_data = cell.getData( rhs.getCellDataID() )->getPointer( level );
            auto dst_data = cell.getData( dst.getCellDataID() )->getPointer( level );
            auto stencil  = cell.getData( cellStencilID_ )->getData( level );

            if ( backwards )
            {
               vertexdof::macrocell::generated::sor_3D_macrocell_P1_backwards(
                   dst_data, rhs_data, static_cast< int32_t >( level ), stencil, relax );
            }
            else
            {
               vertexdof::macrocell::generated::sor_3D_macrocell_P1(
                   dst_data, rhs_data, static_cast< int32_t >( level ), stencil, relax );
            }
         }
         else
         {
            vertexdof::macrocell::smooth_sor< real_t >(
                level, cell, cellStencilID_, dst.getCellDataID(), rhs.getCellDataID(), relax );
         }
      }
   }

   this->timingTree_->stop( "Macro-Cell" );
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::smooth_sor( const P1Function< real_t >& dst,
                                                                                 const P1Function< real_t >& rhs,
                                                                                 real_t                      relax,
                                                                                 size_t                      level,
                                                                                 DoFType                     flag,
                                                                                 const bool&                 backwards ) const
{
   if ( backwards )
   {
      WALBERLA_CHECK( globalDefines::useGeneratedKernels, "Backward SOR only implemented in generated kernels." )
      this->startTiming( "SOR backwards" );
   }
   else
   {
      this->startTiming( "SOR" );
   }

   dst.communicate< Vertex, Edge >( level );
   dst.communicate< Edge, Face >( level );
   dst.communicate< Face, Cell >( level );

   dst.communicate< Cell, Face >( level );
   dst.communicate< Face, Edge >( level );
   dst.communicate< Edge, Vertex >( level );

   if ( backwards )
   {
      smooth_sor_macro_cells( dst, rhs, relax, level, flag, backwards );

      dst.communicate< Cell, Face >( level );

      smooth_sor_macro_faces( dst, rhs, relax, level, flag, backwards );

      dst.communicate< Face, Edge >( level );

      smooth_sor_macro_edges( dst, rhs, relax, level, flag, backwards );

      dst.communicate< Edge, Vertex >( level );

      smooth_sor_macro_vertices( dst, rhs, relax, level, flag, backwards );
   }
   else
   {
      smooth_sor_macro_vertices( dst, rhs, relax, level, flag, backwards );

      dst.communicate< Vertex, Edge >( level );

      smooth_sor_macro_edges( dst, rhs, relax, level, flag, backwards );

      dst.communicate< Edge, Face >( level );

      smooth_sor_macro_faces( dst, rhs, relax, level, flag, backwards );

      dst.communicate< Face, Cell >( level );

      smooth_sor_macro_cells( dst, rhs, relax, level, flag, backwards );
   }

   if ( backwards )
      this->stopTiming( "SOR backwards" );
   else
      this->stopTiming( "SOR" );
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::smooth_jac( const P1Function< real_t >& dst,
                                                                                 const P1Function< real_t >& rhs,
                                                                                 const P1Function< real_t >& src,
                                                                                 real_t                      omega,
                                                                                 size_t                      level,
                                                                                 DoFType                     flag ) const
{
   this->startTiming( "smooth_jac" );

   // compute the current residual
   this->apply( src, dst, level, flag );
   dst.assign( { real_c( 1 ), real_c( -1 ) }, { rhs, dst }, level, flag );

   // perform Jacobi update step
   dst.multElementwise( { *getInverseDiagonalValues(), dst }, level, flag );
   dst.assign( { 1.0, omega }, { src, dst }, level, flag );

   this->stopTiming( "smooth_jac" );
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::scale( real_t scalar )
{
   WALBERLA_CHECK_GREATER_EQUAL( minLevel_, 2, "scale() not implemented for level < 2" )
   WALBERLA_CHECK( !this->storage_->hasGlobalCells(), "scale() not implemented for macro-cells" )

   for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
   {
      for ( auto& it : storage_->getFaces() )
      {
         Face& face         = *it.second;
         auto  face_stencil = face.getData( faceStencilID_ )->getPointer( level );

         for ( const auto& neighbor : vertexdof::macroface::neighborsWithCenter )
         {
            face_stencil[vertexdof::stencilIndexFromVertex( neighbor )] *= scalar;
         }
      }

      for ( auto& it : storage_->getEdges() )
      {
         Edge& edge         = *it.second;
         auto  edge_stencil = edge.getData( edgeStencilID_ )->getPointer( level );

         edge_stencil[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] *= scalar;

         for ( const auto& neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
         {
            edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] *= scalar;
         }

         for ( const auto& neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
         {
            edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] *= scalar;
         }

         if ( edge.getNumNeighborFaces() == 2 )
         {
            for ( const auto& neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
            {
               edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] *= scalar;
            }
         }
      }

      for ( auto& it : storage_->getVertices() )
      {
         Vertex& vertex         = *it.second;
         auto    vertex_stencil = vertex.getData( vertexStencilID_ )->getPointer( level );
         for ( uint_t i = 0; i < vertex.getData( vertexStencilID_ )->getSize( level ); ++i )
         {
            vertex_stencil[i] *= scalar;
         }
      }
   }
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::computeDiagonalOperatorValues( bool invert )
{
   std::shared_ptr< P1Function< real_t > > targetFunction;
   if ( invert )
   {
      if ( !inverseDiagonalValues_ )
      {
         inverseDiagonalValues_ =
             std::make_shared< P1Function< real_t > >( "inverse diagonal entries", storage_, minLevel_, maxLevel_ );
      }
      targetFunction = inverseDiagonalValues_;
   }
   else
   {
      if ( !diagonalValues_ )
      {
         diagonalValues_ = std::make_shared< P1Function< real_t > >( "diagonal entries", storage_, minLevel_, maxLevel_ );
      }
      targetFunction = diagonalValues_;
   }

   for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      for ( const auto& it : storage_->getVertices() )
      {
         auto vertex        = it.second;
         auto stencilMemory = vertex->getData( getVertexStencilID() )->getPointer( level );
         auto targetMemory  = vertex->getData( targetFunction->getVertexDataID() )->getPointer( level );
         targetMemory[0]    = stencilMemory[0];
      }

      if ( level >= 1 )
      {
         for ( const auto& it : storage_->getEdges() )
         {
            auto edge          = it.second;
            auto stencilMemory = edge->getData( getEdgeStencilID() )->getPointer( level );
            auto stencilMap    = edge->getData( getEdgeStencil3DID() )->getData( level );
            auto targetMemory  = edge->getData( targetFunction->getEdgeDataID() )->getPointer( level );

            real_t centerValue = 0;
            if ( storage_->hasGlobalCells() )
            {
               for ( uint_t neighborCellID = 0; neighborCellID < edge->getNumNeighborCells(); neighborCellID++ )
               {
                  for ( auto stencilIt : stencilMap[neighborCellID] )
                  {
                     if ( stencilIt.first == indexing::IndexIncrement( { 0, 0, 0 } ) )
                     {
                        centerValue += stencilIt.second;
                     }
                  }
               }
            }
            else
            {
               centerValue = stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )];
            }

            for ( auto idx : vertexdof::macroedge::Iterator( level ) )
            {
               targetMemory[vertexdof::macroedge::index( level, idx.x() )] = centerValue;
            }
         }
      }

      if ( level >= 1 )
      {
         for ( const auto& it : storage_->getFaces() )
         {
            auto face          = it.second;
            auto stencilMemory = face->getData( getFaceStencilID() )->getPointer( level );
            auto stencilMap    = face->getData( getFaceStencil3DID() )->getData( level );
            auto targetMemory  = face->getData( targetFunction->getFaceDataID() )->getPointer( level );

            real_t centerValue = 0;
            if ( storage_->hasGlobalCells() )
            {
               for ( uint_t neighborCellID = 0; neighborCellID < face->getNumNeighborCells(); neighborCellID++ )
               {
                  for ( auto stencilIt : stencilMap[neighborCellID] )
                  {
                     if ( stencilIt.first == indexing::IndexIncrement( { 0, 0, 0 } ) )
                     {
                        centerValue += stencilIt.second;
                     }
                  }
               }
            }
            else
            {
               centerValue = stencilMemory[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];
            }

            for ( auto idx : vertexdof::macroface::Iterator( level ) )
            {
               targetMemory[vertexdof::macroface::index( level, idx.x(), idx.y() )] = centerValue;
            }
         }
      }

      if ( level >= 2 )
      {
         for ( const auto& it : storage_->getCells() )
         {
            auto cell         = it.second;
            auto stencilMap   = cell->getData( getCellStencilID() )->getData( level );
            auto targetMemory = cell->getData( targetFunction->getCellDataID() )->getPointer( level );

            real_t centerValue = stencilMap[indexing::IndexIncrement( { 0, 0, 0 } )];

            for ( auto idx : vertexdof::macrocell::Iterator( level ) )
            {
               targetMemory[vertexdof::macrocell::index( level, idx.x(), idx.y(), idx.z() )] = centerValue;
            }
         }
      }

      if ( invert )
      {
         targetFunction->invertElementwise( level, All, false );
      }
   }
}

template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, fenics::NoAssemble > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, fenics::UndefinedAssembly > >;

template class P1ConstantOperator<
    P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, fenics::UndefinedAssembly >, true >;

template class P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_2_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_3_otherwise > >;

template class P1ConstantOperator< P1FenicsForm< p1_div_cell_integral_0_otherwise, p1_tet_div_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< p1_div_cell_integral_1_otherwise, p1_tet_div_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_div_tet_cell_integral_2_otherwise > >;

template class P1ConstantOperator< P1FenicsForm< p1_divt_cell_integral_0_otherwise, p1_tet_divt_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< p1_divt_cell_integral_1_otherwise, p1_tet_divt_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_divt_tet_cell_integral_2_otherwise > >;

template class P1ConstantOperator< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise >,
                                   false,
                                   true,
                                   false >;
template class P1ConstantOperator< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise >,
                                   false,
                                   true,
                                   true >;

template class P1ConstantOperator< P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise >,
                                   true,
                                   false,
                                   true >;

template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::UndefinedAssembly, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >;

template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::UndefinedAssembly, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > >;

template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_tet_diffusion_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_tet_mass_cell_integral_0_otherwise > >;

template class P1ConstantOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >;
template class P1ConstantOperator<
    P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >;

template class P1ConstantOperator< P2FenicsForm< p2_divt_cell_integral_0_otherwise, p2_tet_divt_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< p2_divt_cell_integral_1_otherwise, p2_tet_divt_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< fenics::NoAssemble, p2_tet_divt_tet_cell_integral_2_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< p2_div_cell_integral_0_otherwise, p2_tet_div_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< p2_div_cell_integral_1_otherwise, p2_tet_div_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< fenics::NoAssemble, p2_tet_div_tet_cell_integral_2_otherwise > >;

template class P1ConstantOperator<
    P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_0_otherwise, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator<
    P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_1_otherwise, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P1ToP2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > >;

template class P1ConstantOperator<
    P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator<
    P2ToP1FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P2ToP1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >;

template class P1ConstantOperator< P2FenicsForm< p2_pspg_cell_integral_0_otherwise, p2_tet_pspg_tet_cell_integral_0_otherwise > >;

template class P1ConstantOperator< P1LinearCombinationForm >;
template class P1ConstantOperator< P2LinearCombinationForm >;
template class P1ConstantOperator< P2RowSumForm >;

// Mostly for testing, as the P1ConstantOperator inherently was designed to support lumping
template class P1ConstantOperator< P1RowSumForm >;

// The following instantiations are required as building blocks in the P2ConstantEpsilon operator class
// clang-format off
template class P1ConstantOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_0_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_1_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_2_otherwise > >;

template class P1ConstantOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_2_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_3_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_3_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_4_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_5_otherwise > >;

template class P1ConstantOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_6_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_7_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_8_otherwise > >;
// clang-format on

// The following instantiations are required as building blocks in the P2ConstantFullViscousOperator class
// clang-format off
template class P1ConstantOperator< P2FenicsForm< p2_stokes_full_cell_integral_0_otherwise, p2_tet_stokes_full_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< p2_stokes_full_cell_integral_1_otherwise, p2_tet_stokes_full_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_2_otherwise > >;

template class P1ConstantOperator< P2FenicsForm< p2_stokes_full_cell_integral_2_otherwise, p2_tet_stokes_full_tet_cell_integral_3_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< p2_stokes_full_cell_integral_3_otherwise, p2_tet_stokes_full_tet_cell_integral_4_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_5_otherwise > >;

template class P1ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_6_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_7_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_8_otherwise > >;
// clang-format on

} // namespace hyteg
