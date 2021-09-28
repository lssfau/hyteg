/*
 * Copyright (c) 2017-2021 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "EdgeDoFToVertexDoFOperator.hpp"

#include "core/OpenMP.h"

#include "hyteg/forms/P2LinearCombinationForm.hpp"
#include "hyteg/forms/P2RowSumForm.hpp"
#include "hyteg/forms/form_fenics_base/P2ToP1FenicsForm.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFPetsc.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/generatedKernels/apply_2D_macroface_edgedof_to_vertexdof_add.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/generatedKernels/apply_2D_macroface_edgedof_to_vertexdof_replace.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/generatedKernels/apply_3D_macrocell_edgedof_to_vertexdof_add.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/generatedKernels/apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/generatedKernels/apply_3D_macroface_one_sided_edgedof_to_vertexdof_add.hpp"
#include "hyteg/p2functionspace/P2Elements.hpp"
#include "hyteg/p2functionspace/variablestencil/P2VariableStencilCommon.hpp"

#include "EdgeDoFToVertexDoFApply.hpp"

namespace hyteg {

using walberla::int_c;

template < class EdgeDoFToVertexDoFForm >
EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::EdgeDoFToVertexDoFOperator(
    const std::shared_ptr< PrimitiveStorage >& storage,
    const size_t&                              minLevel,
    const size_t&                              maxLevel )
: EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >( storage, minLevel, maxLevel, EdgeDoFToVertexDoFForm() )
{}

template < class EdgeDoFToVertexDoFForm >
EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::EdgeDoFToVertexDoFOperator(
    const std::shared_ptr< PrimitiveStorage >& storage,
    const size_t&                              minLevel,
    const size_t&                              maxLevel,
    const EdgeDoFToVertexDoFForm&              form )
: Operator( storage, minLevel, maxLevel )
, form_( form )
{
   using namespace EdgeDoFToVertexDoF;

   auto vertexDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Vertex > >(
       minLevel_, maxLevel_, macroVertexEdgeDoFToVertexDoFStencilSize );

   auto vertex3DDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< EdgeDoFToVertexDoF::MacroVertexStencilMap_T >, Vertex > >(
           minLevel_, maxLevel_ );

   auto edgeDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Edge > >(
       minLevel_, maxLevel_, macroEdgeEdgeDoFToVertexDoFStencilSize );

   auto edge3DDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge > >(
           minLevel_, maxLevel_ );

   auto faceDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Face > >(
       minLevel_, maxLevel_, macroFaceEdgeDoFToVertexDoFStencilSize );

   auto face3DDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< EdgeDoFToVertexDoF::MacroFaceStencilMap_T >, Face > >(
           minLevel_, maxLevel_ );

   auto cellDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell > >(
           minLevel_, maxLevel_ );

   storage->addVertexData( vertexStencilID_, vertexDataHandling, "VertexDoFToEdgeDoFOperatorVertexStencil" );
   storage->addVertexData( vertexStencil3DID_, vertex3DDataHandling, "VertexDoFToEdgeDoFOperatorVertexStencil3D" );
   storage->addEdgeData( edgeStencilID_, edgeDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil" );
   storage->addEdgeData( edgeStencil3DID_, edge3DDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil3D" );
   storage->addFaceData( faceStencilID_, faceDataHandling, "VertexDoFToEdgeDoFOperatorFaceStencil" );
   storage->addFaceData( faceStencil3DID_, face3DDataHandling, "VertexDoFToEdgeDoFOperatorFaceStencil3D" );
   storage->addCellData( cellStencilID_, cellDataHandling, "VertexDoFToEdgeDoFOperatorCellStencil" );

   if ( this->getStorage()->hasGlobalCells() )
   {
      if ( form_.assemble3D() )
      {
         // WALBERLA_ABORT( "assembleEdgeToVertexStencils< UFCOperator3D > not implemented!" );
         assembleEdgeToVertexStencils< EdgeDoFToVertexDoFForm >( this->getStorage(),
                                                                 this->minLevel_,
                                                                 this->maxLevel_,
                                                                 getVertexStencil3DID(),
                                                                 getEdgeStencil3DID(),
                                                                 getFaceStencil3DID(),
                                                                 getCellStencilID(),
                                                                 form_ );
      }
   }
   else
   {
      // Only assemble stencils if UFCOperator is specified
      if ( form_.assemble2D() )
      {
         assembleStencils();
      }
   }
}

template < class EdgeDoFToVertexDoFForm >
void EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                                     const EdgeDoFFunction< matIdx_t >&          src,
                                                                     const P1Function< matIdx_t >&               dst,
                                                                     size_t                                      level,
                                                                     DoFType                                     flag ) const
{
   auto storage = this->getStorage();

   for ( auto& it : this->getStorage()->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         if ( storage->hasGlobalCells() )
         {
            EdgeDoFToVertexDoF::saveVertexOperator3D(
                level, vertex, *storage, this->getVertexStencil3DID(), src.getVertexDataID(), dst.getVertexDataID(), mat );
         }
         else
         {
            EdgeDoFToVertexDoF::saveVertexOperator(
                level, vertex, this->getVertexStencilID(), src.getVertexDataID(), dst.getVertexDataID(), mat );
         }
      }
   }

   if ( level >= 1 )
   {
      for ( auto& it : this->getStorage()->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            if ( storage->hasGlobalCells() )
            {
               EdgeDoFToVertexDoF::saveEdgeOperator3D(
                   level, edge, *storage, this->getEdgeStencil3DID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat );
            }
            else
            {
               EdgeDoFToVertexDoF::saveEdgeOperator(
                   level, edge, this->getEdgeStencilID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat );
            }
         }
      }
   }

   if ( level >= 2 )
   {
      for ( auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            if ( storage->hasGlobalCells() )
            {
               EdgeDoFToVertexDoF::saveFaceOperator3D(
                   level, face, *storage, this->getFaceStencil3DID(), src.getFaceDataID(), dst.getFaceDataID(), mat );
            }
            else
            {
               EdgeDoFToVertexDoF::saveFaceOperator(
                   level, face, this->getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat );
            }
         }
      }

      for ( auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
         if ( testFlag( cellBC, flag ) )
         {
            EdgeDoFToVertexDoF::saveCellOperator(
                level, cell, this->getCellStencilID(), src.getCellDataID(), dst.getCellDataID(), mat );
         }
      }
   }
}

template < class EdgeDoFToVertexDoFForm >
void EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::assembleStencils()
{
   using namespace P2Elements;

   // Assemble stencils on all levels
   for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
   {
      // Assemble face stencils
      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         uint_t rowsize = levelinfo::num_microvertices_per_edge( level );

         Point3D x( face.coords[0] );
         real_t  h = 1.0 / ( walberla::real_c( rowsize - 1 ) );

         Point3D d0 = h * ( face.coords[1] - face.coords[0] );
         Point3D d2 = h * ( face.coords[2] - face.coords[0] );

         form_.setGeometryMap( face.getGeometryMap() );

         Point3D dirS  = -1.0 * d2;
         Point3D dirSE = d0 - 1.0 * d2;
         Point3D dirE  = d0;
         Point3D dirW  = -1.0 * d0;
         Point3D dirNW = -1.0 * d0 + d2;
         Point3D dirN  = d2;

         real_t* vStencil = storage_->getFace( face.getID() )->getData( faceStencilID_ )->getPointer( level );

         P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
             form_, { x, x + dirW, x + dirS }, P2Elements::P2Face::elementSW_reord, vStencil );
         P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
             form_, { x, x + dirS, x + dirSE }, P2Elements::P2Face::elementS_reord, vStencil );
         P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
             form_, { x, x + dirSE, x + dirE }, P2Elements::P2Face::elementSE_reord, vStencil );
         P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
             form_, { x, x + dirE, x + dirN }, P2Elements::P2Face::elementNE_reord, vStencil );
         P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
             form_, { x, x + dirN, x + dirNW }, P2Elements::P2Face::elementN_reord, vStencil );
         P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
             form_, { x, x + dirNW, x + dirW }, P2Elements::P2Face::elementNW_reord, vStencil );
      }

      // Assemble edge stencils
      for ( auto& it : storage_->getEdges() )
      {
         Edge&   edge     = *it.second;
         real_t* vStencil = storage_->getEdge( edge.getID() )->getData( edgeStencilID_ )->getPointer( level );

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

         Point3D x  = edge.getCoordinates()[0];
         Point3D dx = h * edge.getDirection();
         x += dx;

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
         P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
             form_, { x, x + dir_W, x + dir_S }, P2Elements::P2Face::elementSW_reord, vStencil );
         P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
             form_, { x, x + dir_S, x + dir_SE }, P2Elements::P2Face::elementS_reord, vStencil );
         P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
             form_, { x, x + dir_SE, x + dir_E }, P2Elements::P2Face::elementSE_reord, vStencil );

         if ( edge.getNumNeighborFaces() == 2 )
         {
            form_.setGeometryMap( faceN->getGeometryMap() );
            P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
                form_, { x, x + dir_E, x + dir_N }, P2Elements::P2Face::elementNE_reord, vStencil );
            P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
                form_, { x, x + dir_N, x + dir_NW }, P2Elements::P2Face::elementN_reord, vStencil );
            P2::variablestencil::assembleEdgeToVertexStencil< EdgeDoFToVertexDoFForm >(
                form_, { x, x + dir_NW, x + dir_W }, P2Elements::P2Face::elementNW_reord, vStencil );
         }
      }

      for ( auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         // Assemble EdgeToVertex
         real_t* vStencil = storage_->getVertex( vertex.getID() )->getData( vertexStencilID_ )->getPointer( level );

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
            form_.integrateEdgeToVertex( { { x, x + d0, x + d2 } }, matrixRow );

            uint_t i = 1;
            // iterate over adjacent edges
            for ( auto& edgeId : adj_edges )
            {
               uint_t edge_idx = vertex.edge_index( edgeId );
               vStencil[edge_idx] += matrixRow[3 - i];
               i += 1;
            }

            walberla::uint_t face_idx = vertex.getNumNeighborEdges() + vertex.face_index( face->getID() );
            vStencil[face_idx] += matrixRow[0];

            ++neighborId;
         }
      }
   }
}

template < class EdgeDoFToVertexDoFForm >
void EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::apply( const EdgeDoFFunction< real_t >& src,
                                                                  const P1Function< double >&      dst,
                                                                  uint_t                           level,
                                                                  DoFType                          flag,
                                                                  UpdateType                       updateType ) const
{
   using namespace EdgeDoFToVertexDoF;
   this->startTiming( "Apply" );

   ///there might be room for optimization in the communication. i.e. splitting communicate into start and end to overlap comm and calc

   src.communicate< Face, Cell >( level );

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
            if ( hyteg::globalDefines::useGeneratedKernels && updateType == Add )
            {
               typedef edgedof::EdgeDoFOrientation eo;
               auto                                dstData     = cell.getData( dst.getCellDataID() )->getPointer( level );
               auto                                srcData     = cell.getData( src.getCellDataID() )->getPointer( level );
               auto                                stencilData = cell.getData( cellStencilID_ )->getData( level );
               std::map< eo, uint_t >              firstIdx;
               for ( auto e : edgedof::allEdgeDoFOrientations )
                  firstIdx[e] = edgedof::macrocell::index( level, 0, 0, 0, e );
               EdgeDoFToVertexDoF::generated::apply_3D_macrocell_edgedof_to_vertexdof_add( &srcData[firstIdx[eo::X]],
                                                                                           &srcData[firstIdx[eo::XY]],
                                                                                           &srcData[firstIdx[eo::XYZ]],
                                                                                           &srcData[firstIdx[eo::XZ]],
                                                                                           &srcData[firstIdx[eo::Y]],
                                                                                           &srcData[firstIdx[eo::YZ]],
                                                                                           &srcData[firstIdx[eo::Z]],
                                                                                           dstData,
                                                                                           stencilData,
                                                                                           static_cast< int32_t >( level ) );
            }
            else
            {
               applyCell( level, cell, cellStencilID_, src.getCellDataID(), dst.getCellDataID(), updateType );
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Cell" );

   src.communicate< Edge, Face >( level );
   src.communicate< Cell, Face >( level );

   this->timingTree_->start( "Macro-Face" );

   if ( level >= 1 )
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
               if ( hyteg::globalDefines::useGeneratedKernels && updateType == Add )
               {
                  WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start( "Generated" ); }
                  auto dstData     = face.getData( dst.getFaceDataID() )->getPointer( level );
                  auto srcData     = face.getData( src.getFaceDataID() )->getPointer( level );
                  auto stencilData = face.getData( faceStencil3DID_ )->getData( level );

                  typedef edgedof::EdgeDoFOrientation eo;

                  std::map< eo, uint_t > firstIdxInner;
                  for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
                  {
                     firstIdxInner[e] = edgedof::macroface::index( level, 0, 0, e );
                  }

                  std::map< uint_t, std::map< eo, uint_t > > firstIdxNeighbor;
                  for ( uint_t neighbor = 0; neighbor < 2; neighbor++ )
                  {
                     for ( auto e : edgedof::allEdgeDoFOrientations )
                     {
                        firstIdxNeighbor[neighbor][e] = edgedof::macroface::index( level, 0, 0, e, neighbor );
                     }
                  }

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

                  EdgeDoFToVertexDoF::generated::apply_3D_macroface_one_sided_edgedof_to_vertexdof_add(
                      &srcData[firstIdxInner[eo::X]],
                      &srcData[firstIdxInner[eo::XY]],
                      &srcData[firstIdxInner[eo::Y]],
                      &srcData[firstIdxNeighbor[0][eo::X]],
                      &srcData[firstIdxNeighbor[0][eo::XY]],
                      &srcData[firstIdxNeighbor[0][eo::XYZ]],
                      &srcData[firstIdxNeighbor[0][eo::XZ]],
                      &srcData[firstIdxNeighbor[0][eo::Y]],
                      &srcData[firstIdxNeighbor[0][eo::YZ]],
                      &srcData[firstIdxNeighbor[0][eo::Z]],
                      dstData,
                      stencilData[0],
                      static_cast< int32_t >( level ),
                      neighbor_cell_0_local_vertex_id_0,
                      neighbor_cell_0_local_vertex_id_1,
                      neighbor_cell_0_local_vertex_id_2 );

                  if ( face.getNumNeighborCells() == 2 )
                  {
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

                     EdgeDoFToVertexDoF::generated::apply_3D_macroface_one_sided_edgedof_to_vertexdof_add(
                         &srcData[firstIdxInner[eo::X]],
                         &srcData[firstIdxInner[eo::XY]],
                         &srcData[firstIdxInner[eo::Y]],
                         &srcData[firstIdxNeighbor[1][eo::X]],
                         &srcData[firstIdxNeighbor[1][eo::XY]],
                         &srcData[firstIdxNeighbor[1][eo::XYZ]],
                         &srcData[firstIdxNeighbor[1][eo::XZ]],
                         &srcData[firstIdxNeighbor[1][eo::Y]],
                         &srcData[firstIdxNeighbor[1][eo::YZ]],
                         &srcData[firstIdxNeighbor[1][eo::Z]],
                         dstData,
                         stencilData[1],
                         static_cast< int32_t >( level ),
                         neighbor_cell_1_local_vertex_id_0,
                         neighbor_cell_1_local_vertex_id_1,
                         neighbor_cell_1_local_vertex_id_2 );
                  }
                  WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->stop( "Generated" ); }
               }
               else
               {
                  WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start( "Not generated" ); }
                  applyFace3D( level, face, *storage_, faceStencil3DID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
                  WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->stop( "Not generated" ); }
               }
            }
            else
            {
               if ( hyteg::globalDefines::useGeneratedKernels )
               {
                  real_t* opr_data = face.getData( faceStencilID_ )->getPointer( level );
                  real_t* src_data = face.getData( src.getFaceDataID() )->getPointer( level );
                  real_t* dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );

                  typedef edgedof::EdgeDoFOrientation eo;
                  std::map< eo, uint_t >              firstIdx;
                  for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
                     firstIdx[e] = edgedof::macroface::index( level, 0, 0, e );

                  if ( updateType == hyteg::Replace )
                  {
                     EdgeDoFToVertexDoF::generated::apply_2D_macroface_edgedof_to_vertexdof_replace(
                         &src_data[firstIdx[eo::X]],
                         &src_data[firstIdx[eo::XY]],
                         &src_data[firstIdx[eo::Y]],
                         opr_data,
                         dst_data,
                         static_cast< int32_t >( level ) );
                  }
                  else if ( updateType == hyteg::Add )
                  {
                     EdgeDoFToVertexDoF::generated::apply_2D_macroface_edgedof_to_vertexdof_add(
                         &src_data[firstIdx[eo::X]],
                         &src_data[firstIdx[eo::XY]],
                         &src_data[firstIdx[eo::Y]],
                         opr_data,
                         dst_data,
                         static_cast< int32_t >( level ) );
                  }
               }
               else
               {
                  applyFace( level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
               }
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Face" );

   src.communicate< Face, Edge >( level );

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
            if ( storage_->hasGlobalCells() )
            {
               if ( globalDefines::useGeneratedKernels && updateType == Add )
               {
                  auto dstData     = edge.getData( dst.getEdgeDataID() )->getPointer( level );
                  auto srcData     = edge.getData( src.getEdgeDataID() )->getPointer( level );
                  auto stencilData = edge.getData( edgeStencil3DID_ )->getData( level );

                  for ( uint_t cellID = 0; cellID < edge.getNumNeighborCells(); cellID++ )
                  {
                     const auto neighborCell    = storage_->getCell( edge.neighborCells().at( cellID ) );
                     const auto cellLocalEdgeID = neighborCell->getLocalEdgeID( edge.getID() );

                     const auto cellLocalVertexID0 =
                         neighborCell->getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 0 );
                     const auto cellLocalVertexID1 =
                         neighborCell->getEdgeLocalVertexToCellLocalVertexMaps().at( cellLocalEdgeID ).at( 1 );
                     const auto cellLocalVertexID2 =
                         algorithms::getMissingIntegersAscending< 2, 4 >( { { cellLocalVertexID0, cellLocalVertexID1 } } )
                             .at( 2 );

                     const std::vector< uint_t > neighborFaces(
                         indexing::cellLocalEdgeIDsToCellLocalNeighborFaceIDs.at( cellLocalEdgeID ).begin(),
                         indexing::cellLocalEdgeIDsToCellLocalNeighborFaceIDs.at( cellLocalEdgeID ).end() );
                     const auto cellLocalNeighborFaceID0 =
                         neighborFaces.at( 0 ) < neighborFaces.at( 1 ) ? neighborFaces.at( 0 ) : neighborFaces.at( 1 );
                     const auto cellLocalNeighborFaceID1 =
                         neighborFaces.at( 0 ) > neighborFaces.at( 1 ) ? neighborFaces.at( 0 ) : neighborFaces.at( 1 );

                     const auto neighborFacePrimitiveID0 = neighborCell->neighborFaces().at( cellLocalNeighborFaceID0 );
                     const auto neighborFacePrimitiveID1 = neighborCell->neighborFaces().at( cellLocalNeighborFaceID1 );

                     const auto edgeLocalFaceID0 = edge.face_index( neighborFacePrimitiveID0 );
                     const auto edgeLocalFaceID1 = edge.face_index( neighborFacePrimitiveID1 );

                     EdgeDoFToVertexDoF::generated::apply_3D_macroedge_per_cell_edgedof_to_vertexdof_add(
                         srcData,
                         dstData,
                         static_cast< int64_t >( cellID ),
                         stencilData[uint_c( cellID )],
                         static_cast< int64_t >( edgeLocalFaceID0 ),
                         static_cast< int64_t >( edgeLocalFaceID1 ),
                         static_cast< int32_t >( level ),
                         static_cast< int64_t >( cellLocalVertexID0 ),
                         static_cast< int64_t >( cellLocalVertexID1 ),
                         static_cast< int64_t >( cellLocalVertexID2 ),
                         static_cast< int64_t >( edge.getNumNeighborFaces() ) );
                  }
               }
               else
               {
                  applyEdge3D(
                      level, edge, *getStorage(), edgeStencil3DID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
               }
            }
            else
            {
               applyEdge( level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Edge" );

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
         if ( storage_->hasGlobalCells() )
         {
            applyVertex3D(
                level, vertex, *getStorage(), vertexStencil3DID_, src.getVertexDataID(), dst.getVertexDataID(), updateType );
         }
         else
         {
            applyVertex( level, vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), updateType );
         }
      }
   }

   this->timingTree_->stop( "Macro-Vertex" );

   this->stopTiming( "Apply" );
}

template < class EdgeDoFToVertexDoFForm >
const PrimitiveDataID< StencilMemory< real_t >, Vertex >&
    EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::getVertexStencilID() const
{
   return vertexStencilID_;
}

template < class EdgeDoFToVertexDoFForm >
const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroVertexStencilMap_T >, Vertex >&
    EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::getVertexStencil3DID() const
{
   return vertexStencil3DID_;
}

template < class EdgeDoFToVertexDoFForm >
const PrimitiveDataID< StencilMemory< real_t >, Edge >&
    EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::getEdgeStencilID() const
{
   return edgeStencilID_;
}

template < class EdgeDoFToVertexDoFForm >
const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroEdgeStencilMap_T >, Edge >&
    EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::getEdgeStencil3DID() const
{
   return edgeStencil3DID_;
}

template < class EdgeDoFToVertexDoFForm >
const PrimitiveDataID< StencilMemory< real_t >, Face >&
    EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::getFaceStencilID() const
{
   return faceStencilID_;
}

template < class EdgeDoFToVertexDoFForm >
const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroFaceStencilMap_T >, Face >&
    EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::getFaceStencil3DID() const
{
   return faceStencil3DID_;
}

template < class EdgeDoFToVertexDoFForm >
const PrimitiveDataID< LevelWiseMemory< EdgeDoFToVertexDoF::MacroCellStencilMap_T >, Cell >&
    EdgeDoFToVertexDoFOperator< EdgeDoFToVertexDoFForm >::getCellStencilID() const
{
   return cellStencilID_;
}

namespace EdgeDoFToVertexDoF {
////////// Stencil sizes //////////
uint_t macroVertexEdgeDoFToVertexDoFStencilSize( const uint_t& level, const Primitive& primitive )
{
   WALBERLA_UNUSED( level );
   return primitive.getNumNeighborEdges() + primitive.getNumNeighborFaces();
}

uint_t macroEdgeEdgeDoFToVertexDoFStencilSize( const uint_t& level, const Primitive& primitive )
{
   WALBERLA_UNUSED( level );
   return 2 + 5 * primitive.getNumNeighborFaces();
}

uint_t macroFaceEdgeDoFToVertexDoFStencilSize( const uint_t& level, const Primitive& primitive )
{
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( primitive );
   return 12;
}

uint_t macroCellEdgeDoFToVertexDoFStencilSize( const uint_t& level, const Primitive& primitive )
{
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( primitive );
   return 7 * 27;
}

} // namespace EdgeDoFToVertexDoF

template class EdgeDoFToVertexDoFOperator< P2FenicsForm< hyteg::fenics::NoAssemble, hyteg::fenics::NoAssemble > >;
template class EdgeDoFToVertexDoFOperator<
    P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >;
template class EdgeDoFToVertexDoFOperator<
    P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >;

template class EdgeDoFToVertexDoFOperator<
    P2FenicsForm< p2_divt_cell_integral_0_otherwise, p2_tet_divt_tet_cell_integral_0_otherwise > >;
template class EdgeDoFToVertexDoFOperator<
    P2FenicsForm< p2_divt_cell_integral_1_otherwise, p2_tet_divt_tet_cell_integral_1_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< fenics::NoAssemble, p2_tet_divt_tet_cell_integral_2_otherwise > >;
template class EdgeDoFToVertexDoFOperator<
    P2FenicsForm< p2_div_cell_integral_0_otherwise, p2_tet_div_tet_cell_integral_0_otherwise > >;
template class EdgeDoFToVertexDoFOperator<
    P2FenicsForm< p2_div_cell_integral_1_otherwise, p2_tet_div_tet_cell_integral_1_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< fenics::NoAssemble, p2_tet_div_tet_cell_integral_2_otherwise > >;

template class EdgeDoFToVertexDoFOperator<
    P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > >;
template class EdgeDoFToVertexDoFOperator<
    P2ToP1FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > >;
template class EdgeDoFToVertexDoFOperator<
    P2ToP1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >;

template class EdgeDoFToVertexDoFOperator<
    P2FenicsForm< p2_pspg_cell_integral_0_otherwise, p2_tet_pspg_tet_cell_integral_0_otherwise > >;

template class EdgeDoFToVertexDoFOperator< P2LinearCombinationForm >;
template class EdgeDoFToVertexDoFOperator< P2RowSumForm >;

// The following instantiations are required as building blocks in the P2ConstantEpsilon operator class
// clang-format off
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_0_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_0_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_1_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_1_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_2_otherwise > >;

template class EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_2_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_3_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_3_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_4_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_5_otherwise > >;

template class EdgeDoFToVertexDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_6_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_7_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_8_otherwise > >;
// clang-format on

// The following instantiations are required as building blocks in the P2ConstantFullViscousOperator class
// clang-format off
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_stokes_full_cell_integral_0_otherwise, p2_tet_stokes_full_tet_cell_integral_0_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_stokes_full_cell_integral_1_otherwise, p2_tet_stokes_full_tet_cell_integral_1_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_2_otherwise > >;

template class EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_stokes_full_cell_integral_2_otherwise, p2_tet_stokes_full_tet_cell_integral_3_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< p2_stokes_full_cell_integral_3_otherwise, p2_tet_stokes_full_tet_cell_integral_4_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_5_otherwise > >;

template class EdgeDoFToVertexDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_6_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_7_otherwise > >;
template class EdgeDoFToVertexDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_8_otherwise > >;
// clang-format on

} // namespace hyteg
