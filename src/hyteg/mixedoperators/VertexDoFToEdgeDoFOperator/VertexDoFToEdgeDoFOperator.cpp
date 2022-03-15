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
#include "VertexDoFToEdgeDoFOperator.hpp"

#include "hyteg/forms/P2LinearCombinationForm.hpp"
#include "hyteg/forms/P2RowSumForm.hpp"
#include "hyteg/forms/form_fenics_base/P1ToP2FenicsForm.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFPetsc.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/generatedKernels/apply_2D_macroface_vertexdof_to_edgedof_add.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/generatedKernels/apply_2D_macroface_vertexdof_to_edgedof_replace.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/generatedKernels/apply_3D_macrocell_vertexdof_to_edgedof_add.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/generatedKernels/apply_3D_macrocell_vertexdof_to_edgedof_replace.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/generatedKernels/apply_3D_macroface_one_sided_vertexdof_to_edgedof_add.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/generatedKernels/apply_3D_macroface_one_sided_vertexdof_to_edgedof_replace.hpp"
#include "hyteg/p2functionspace/P2Elements.hpp"
#include "hyteg/p2functionspace/variablestencil/P2VariableStencilCommon.hpp"

namespace hyteg {

using walberla::int_c;

template < class VertexDoFToEdgeDoFForm >
VertexDoFToEdgeDoFOperator< VertexDoFToEdgeDoFForm >::VertexDoFToEdgeDoFOperator(
    const std::shared_ptr< PrimitiveStorage >& storage,
    size_t                                     minLevel,
    size_t                                     maxLevel )
: VertexDoFToEdgeDoFOperator< VertexDoFToEdgeDoFForm >( storage, minLevel, maxLevel, VertexDoFToEdgeDoFForm() )
{}

template < class VertexDoFToEdgeDoFForm >
VertexDoFToEdgeDoFOperator< VertexDoFToEdgeDoFForm >::VertexDoFToEdgeDoFOperator(
    const std::shared_ptr< PrimitiveStorage >& storage,
    size_t                                     minLevel,
    size_t                                     maxLevel,
    const VertexDoFToEdgeDoFForm&              form )
: Operator( storage, minLevel, maxLevel )
, form_( form )
{
   /// since the Vertex does not own any EdgeDoFs only edge and face are needed

   auto edgeDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Edge > >(
       minLevel_, maxLevel_, VertexDoFToEdgeDoF::macroEdgeVertexDoFToEdgeDoFStencilSize );

   auto edge3DDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< VertexDoFToEdgeDoF::MacroEdgeStencilMap_T >, Edge > >(
           minLevel_, maxLevel_ );

   auto faceDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Face > >(
       minLevel_, maxLevel_, VertexDoFToEdgeDoF::macroFaceVertexDoFToEdgeDoFStencilSize );

   auto face3DDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< VertexDoFToEdgeDoF::MacroFaceStencilMap_T >, Face > >(
           minLevel_, maxLevel_ );

   auto cellDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< VertexDoFToEdgeDoF::MacroCellStencilMap_T >, Cell > >(
           minLevel_, maxLevel_ );

   storage->addEdgeData( edgeStencilID_, edgeDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil" );
   storage->addEdgeData( edgeStencil3DID_, edge3DDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil3D" );
   storage->addFaceData( faceStencilID_, faceDataHandling, "VertexDoFToEdgeDoFOperatorFaceStencil" );
   storage->addFaceData( faceStencil3DID_, face3DDataHandling, "VertexDoFToEdgeDoFOperatorFaceStencil3D" );
   storage->addCellData( cellStencilID_, cellDataHandling, "VertexDoFToEdgeDoFOperatorCellStencil" );

   if ( this->getStorage()->hasGlobalCells() )
   {
      if ( form_.assemble3D() )
      {
         // WALBERLA_ABORT("Not implemented.");
         assembleVertexToEdgeStencils< VertexDoFToEdgeDoFForm >( this->getStorage(),
                                                                 this->minLevel_,
                                                                 this->maxLevel_,
                                                                 getEdgeStencil3DID(),
                                                                 getFaceStencil3DID(),
                                                                 getCellStencilID(),
                                                                 form_ );
      }
   }
   else
   {
      if ( form_.assemble2D() )
      {
         assembleStencils();
      }
   }
}

template < class VertexDoFToEdgeDoFForm >
void VertexDoFToEdgeDoFOperator< VertexDoFToEdgeDoFForm >::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                                     const P1Function< idx_t >&                  src,
                                                                     const EdgeDoFFunction< idx_t >&             dst,
                                                                     size_t                                      level,
                                                                     DoFType                                     flag ) const
{
   const auto storage = src.getStorage();

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         if ( storage->hasGlobalCells() )
         {
            VertexDoFToEdgeDoF::saveEdgeOperator3D(
                level, edge, *storage, this->getEdgeStencil3DID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat );
         }
         else
         {
            VertexDoFToEdgeDoF::saveEdgeOperator(
                level, edge, this->getEdgeStencilID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat );
         }
      }
   }

   if ( level >= 1 )
   {
      for ( auto& it : this->getStorage()->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            if ( storage->hasGlobalCells() )
            {
               VertexDoFToEdgeDoF::saveFaceOperator3D(
                   level, face, *storage, this->getFaceStencil3DID(), src.getFaceDataID(), dst.getFaceDataID(), mat );
            }
            else
            {
               VertexDoFToEdgeDoF::saveFaceOperator(
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
            VertexDoFToEdgeDoF::saveCellOperator(
                level, cell, this->getCellStencilID(), src.getCellDataID(), dst.getCellDataID(), mat );
         }
      }
   }
}

template < class VertexDoFToEdgeDoFForm >
void VertexDoFToEdgeDoFOperator< VertexDoFToEdgeDoFForm >::assembleStencils()
{
   using namespace P2Elements;

   // Initialize memory for local 6x6 matrices
   Matrix6r local_stiffness_gray;
   Matrix6r local_stiffness_blue;

   // Assemble stencils on all levels
   for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
   {
      // Assemble face stencils
      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         // Assemble vertexToEdge stencil
         real_t* vStencil = storage_->getFace( face.getID() )->getData( faceStencilID_ )->getPointer( level );

         form_.setGeometryMap( face.getGeometryMap() );

         const Point3D faceBottomLeftCoords  = face.getCoordinates()[0];
         const Point3D faceBottomRightCoords = face.getCoordinates()[1];
         const Point3D faceTopLeftCoords     = face.getCoordinates()[2];

         const Point3D horizontalMicroEdgeOffset = ( ( faceBottomRightCoords - faceBottomLeftCoords ) /
                                                     walberla::real_c( levelinfo::num_microedges_per_edge( level ) ) ) *
                                                   0.5;
         const Point3D verticalMicroEdgeOffset =
             ( ( faceTopLeftCoords - faceBottomLeftCoords ) / walberla::real_c( levelinfo::num_microedges_per_edge( level ) ) ) *
             0.5;

         const Point3D dirHO_W  = -horizontalMicroEdgeOffset;
         const Point3D dirHO_E  = horizontalMicroEdgeOffset;
         const Point3D dirHO_SE = horizontalMicroEdgeOffset - 2.0 * verticalMicroEdgeOffset;
         const Point3D dirHO_NW = -horizontalMicroEdgeOffset + 2.0 * verticalMicroEdgeOffset;

         const Point3D dirVE_N  = verticalMicroEdgeOffset;
         const Point3D dirVE_S  = -verticalMicroEdgeOffset;
         const Point3D dirVE_NW = -2.0 * horizontalMicroEdgeOffset + verticalMicroEdgeOffset;
         const Point3D dirVE_SE = 2.0 * horizontalMicroEdgeOffset - verticalMicroEdgeOffset;

         const Point3D dirDI_SE = horizontalMicroEdgeOffset - verticalMicroEdgeOffset;
         const Point3D dirDI_NE = horizontalMicroEdgeOffset + verticalMicroEdgeOffset;
         const Point3D dirDI_NW = -horizontalMicroEdgeOffset + verticalMicroEdgeOffset;
         const Point3D dirDI_SW = -horizontalMicroEdgeOffset - verticalMicroEdgeOffset;

         auto edgeIt = edgedof::macroface::Iterator( level, 0 );

         // Loop until first interior DoF is reached
         while ( edgeIt->row() == 0 || edgeIt->col() == 0 ||
                 edgeIt->col() + edgeIt->row() == int_c( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
         {
            edgeIt++;
         }

         const Point3D horizontalMicroEdgePosition =
             faceBottomLeftCoords + ( walberla::real_c( edgeIt->col() * 2 + 1 ) * horizontalMicroEdgeOffset +
                                      walberla::real_c( edgeIt->row() * 2 ) * verticalMicroEdgeOffset );
         const Point3D verticalMicroEdgePosition =
             faceBottomLeftCoords + ( walberla::real_c( edgeIt->col() * 2 ) * horizontalMicroEdgeOffset +
                                      walberla::real_c( edgeIt->row() * 2 + 1 ) * verticalMicroEdgeOffset );
         const Point3D diagonalMicroEdgePosition = horizontalMicroEdgePosition + verticalMicroEdgeOffset;

         P2::variablestencil::assembleVertexToEdgeStencil( form_,
                                                           {horizontalMicroEdgePosition + dirHO_W,
                                                            horizontalMicroEdgePosition + dirHO_E,
                                                            horizontalMicroEdgePosition + dirHO_NW},
                                                           {vertexdof::stencilIndexFromHorizontalEdge( SD::VERTEX_W ),
                                                            vertexdof::stencilIndexFromHorizontalEdge( SD::VERTEX_E ),
                                                            vertexdof::stencilIndexFromHorizontalEdge( SD::VERTEX_NW )},
                                                           vStencil );
         P2::variablestencil::assembleVertexToEdgeStencil( form_,
                                                           {horizontalMicroEdgePosition + dirHO_W,
                                                            horizontalMicroEdgePosition + dirHO_E,
                                                            horizontalMicroEdgePosition + dirHO_SE},
                                                           {vertexdof::stencilIndexFromHorizontalEdge( SD::VERTEX_W ),
                                                            vertexdof::stencilIndexFromHorizontalEdge( SD::VERTEX_E ),
                                                            vertexdof::stencilIndexFromHorizontalEdge( SD::VERTEX_SE )},
                                                           vStencil );

         P2::variablestencil::assembleVertexToEdgeStencil(
             form_,
             {verticalMicroEdgePosition + dirVE_N, verticalMicroEdgePosition + dirVE_S, verticalMicroEdgePosition + dirVE_NW},
             {vertexdof::stencilIndexFromVerticalEdge( SD::VERTEX_N ),
              vertexdof::stencilIndexFromVerticalEdge( SD::VERTEX_S ),
              vertexdof::stencilIndexFromVerticalEdge( SD::VERTEX_NW )},
             vStencil );
         P2::variablestencil::assembleVertexToEdgeStencil(
             form_,
             {verticalMicroEdgePosition + dirVE_N, verticalMicroEdgePosition + dirVE_S, verticalMicroEdgePosition + dirVE_SE},
             {vertexdof::stencilIndexFromVerticalEdge( SD::VERTEX_N ),
              vertexdof::stencilIndexFromVerticalEdge( SD::VERTEX_S ),
              vertexdof::stencilIndexFromVerticalEdge( SD::VERTEX_SE )},
             vStencil );

         P2::variablestencil::assembleVertexToEdgeStencil(
             form_,
             {diagonalMicroEdgePosition + dirDI_NW, diagonalMicroEdgePosition + dirDI_SE, diagonalMicroEdgePosition + dirDI_SW},
             {vertexdof::stencilIndexFromDiagonalEdge( SD::VERTEX_NW ),
              vertexdof::stencilIndexFromDiagonalEdge( SD::VERTEX_SE ),
              vertexdof::stencilIndexFromDiagonalEdge( SD::VERTEX_SW )},
             vStencil );
         P2::variablestencil::assembleVertexToEdgeStencil(
             form_,
             {diagonalMicroEdgePosition + dirDI_NW, diagonalMicroEdgePosition + dirDI_SE, diagonalMicroEdgePosition + dirDI_NE},
             {vertexdof::stencilIndexFromDiagonalEdge( SD::VERTEX_NW ),
              vertexdof::stencilIndexFromDiagonalEdge( SD::VERTEX_SE ),
              vertexdof::stencilIndexFromDiagonalEdge( SD::VERTEX_NE )},
             vStencil );
      }

      // Assemble edge stencils
      for ( auto& it : storage_->getEdges() )
      {
         Edge&   edge     = *it.second;
         real_t* vStencil = storage_->getEdge( edge.getID() )->getData( edgeStencilID_ )->getPointer( level );

         using namespace hyteg::edgedof::macroedge;
         size_t rowsize = levelinfo::num_microedges_per_edge( level );

         Face*  faceS   = storage_->getFace( edge.neighborFaces()[0] );
         Face*  faceN   = nullptr;
         uint_t s_south = faceS->vertex_index( edge.neighborVertices()[0] );
         uint_t e_south = faceS->vertex_index( edge.neighborVertices()[1] );
         uint_t o_south = faceS->vertex_index( faceS->get_vertex_opposite_to_edge( edge.getID() ) );

         real_t h = 1.0 / ( walberla::real_c( rowsize ) );

         Point3D dS_se = h * ( faceS->getCoordinates()[e_south] - faceS->getCoordinates()[s_south] );
         //       Point3D dS_so = h * ( faceS->getCoordinates()[o_south] - faceS->getCoordinates()[s_south] );
         Point3D dS_oe = h * ( faceS->getCoordinates()[e_south] - faceS->getCoordinates()[o_south] );

         Point3D dir_SE = 0.5 * dS_se - 1.0 * dS_oe;
         Point3D dir_E  = 0.5 * dS_se;
         Point3D dir_W  = -0.5 * dS_se;

         uint_t  s_north, o_north;
         Point3D dir_NW;

         if ( edge.getNumNeighborFaces() == 2 )
         {
            faceN   = storage_->getFace( edge.neighborFaces()[1] );
            s_north = faceN->vertex_index( edge.neighborVertices()[0] );
            //          e_north = faceN->vertex_index( edge.neighborVertices()[1] );
            o_north = faceN->vertex_index( faceN->get_vertex_opposite_to_edge( edge.getID() ) );

            Point3D dN_so = h * ( faceN->getCoordinates()[o_north] - faceN->getCoordinates()[s_north] );
            //          Point3D dN_oe = h * ( faceN->getCoordinates()[e_north] - faceN->getCoordinates()[o_north] );

            dir_NW = -0.5 * dS_se + 1.0 * dN_so;
         }

         const Point3D leftCoords = edge.getCoordinates()[0];

         std::vector< real_t > vertexToEdge( 4 );
         std::vector< real_t > edgeToEdge( 5 );

         Point3D horizontalMicroEdgePosition;

         horizontalMicroEdgePosition = leftCoords + ( real_c( 0 ) + 0.5 ) * dS_se;

         form_.setGeometryMap( faceS->getGeometryMap() );
         P2::variablestencil::assembleVertexToEdgeStencil(
             form_,
             {horizontalMicroEdgePosition + dir_W, horizontalMicroEdgePosition + dir_E, horizontalMicroEdgePosition + dir_SE},
             {vertexdof::stencilIndexFromHorizontalEdge( SD::VERTEX_W ),
              vertexdof::stencilIndexFromHorizontalEdge( SD::VERTEX_E ),
              vertexdof::stencilIndexFromHorizontalEdge( SD::VERTEX_SE )},
             vStencil );

         if ( edge.getNumNeighborFaces() == 2 )
         {
            form_.setGeometryMap( faceN->getGeometryMap() );
            P2::variablestencil::assembleVertexToEdgeStencil(
                form_,
                {horizontalMicroEdgePosition + dir_W, horizontalMicroEdgePosition + dir_E, horizontalMicroEdgePosition + dir_NW},
                {vertexdof::stencilIndexFromHorizontalEdge( SD::VERTEX_W ),
                 vertexdof::stencilIndexFromHorizontalEdge( SD::VERTEX_E ),
                 vertexdof::stencilIndexFromHorizontalEdge( SD::VERTEX_NW )},
                vStencil );
         }
      }
   }
}

template < class VertexDoFToEdgeDoFForm >
void VertexDoFToEdgeDoFOperator< VertexDoFToEdgeDoFForm >::apply( const P1Function< real_t >&      src,
                                                                  const EdgeDoFFunction< real_t >& dst,
                                                                  size_t                           level,
                                                                  DoFType                          flag,
                                                                  UpdateType                       updateType ) const
{
   this->startTiming( "Apply" );
   ///the order of communication is crucial here.
   ///first the vertex dofs on the macro vertex need to be communicated to the edge since they are needed on the edge and the face
   src.communicate< Vertex, Edge >( level );
   ///secondly the vertex dofs on the macro edge are communicated to the face passing on the vertex dof from the macro vertex
   src.communicate< Edge, Face >( level );
   src.communicate< Face, Cell >( level );
   src.communicate< Cell, Face >( level );
   ///lastly the vertex dofs on the macro face are communicated to the edge which also contain vertex dofs which are located on neighboring edges
   src.startCommunication< Face, Edge >( level );

   this->timingTree_->start( "Macro-Cell" );

   if ( level >= 1 )
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
               typedef edgedof::EdgeDoFOrientation eo;
               auto                                dstData     = cell.getData( dst.getCellDataID() )->getPointer( level );
               auto                                srcData     = cell.getData( src.getCellDataID() )->getPointer( level );
               auto                                stencilData = cell.getData( cellStencilID_ )->getData( level );
               std::map< eo, uint_t >              firstIdx;
               for ( auto e : edgedof::allEdgeDoFOrientations )
                  firstIdx[e] = edgedof::macrocell::index( level, 0, 0, 0, e );

               if ( updateType == Replace )
               {
                  VertexDoFToEdgeDoF::generated::apply_3D_macrocell_vertexdof_to_edgedof_replace( &dstData[firstIdx[eo::X]],
                                                                                                  &dstData[firstIdx[eo::XY]],
                                                                                                  &dstData[firstIdx[eo::XYZ]],
                                                                                                  &dstData[firstIdx[eo::XZ]],
                                                                                                  &dstData[firstIdx[eo::Y]],
                                                                                                  &dstData[firstIdx[eo::YZ]],
                                                                                                  &dstData[firstIdx[eo::Z]],
                                                                                                  srcData,
                                                                                                  static_cast< int32_t >( level ),
                                                                                                  stencilData );
               }
               else
               {
                  VertexDoFToEdgeDoF::generated::apply_3D_macrocell_vertexdof_to_edgedof_add( &dstData[firstIdx[eo::X]],
                                                                                              &dstData[firstIdx[eo::XY]],
                                                                                              &dstData[firstIdx[eo::XYZ]],
                                                                                              &dstData[firstIdx[eo::XZ]],
                                                                                              &dstData[firstIdx[eo::Y]],
                                                                                              &dstData[firstIdx[eo::YZ]],
                                                                                              &dstData[firstIdx[eo::Z]],
                                                                                              srcData,
                                                                                              static_cast< int32_t >( level ),
                                                                                              stencilData );
               }
            }
            else
            {
               VertexDoFToEdgeDoF::applyCell( level, cell, cellStencilID_, src.getCellDataID(), dst.getCellDataID(), updateType );
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Cell" );

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
               if ( hyteg::globalDefines::useGeneratedKernels )
               {
                  auto dstData     = face.getData( dst.getFaceDataID() )->getPointer( level );
                  auto srcData     = face.getData( src.getFaceDataID() )->getPointer( level );
                  auto stencilData = face.getData( faceStencil3DID_ )->getData( level );

                  typedef edgedof::EdgeDoFOrientation eo;
                  std::map< eo, uint_t >              firstIdx;
                  for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
                     firstIdx[e] = edgedof::macroface::index( level, 0, 0, e );

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
                     VertexDoFToEdgeDoF::generated::apply_3D_macroface_one_sided_vertexdof_to_edgedof_replace(
                         &dstData[firstIdx[eo::X]],
                         &dstData[firstIdx[eo::XY]],
                         &dstData[firstIdx[eo::Y]],
                         &srcData[0],
                         &srcData[offset_gl_0],
                         static_cast< int32_t >( level ),
                         neighbor_cell_0_local_vertex_id_0,
                         neighbor_cell_0_local_vertex_id_1,
                         neighbor_cell_0_local_vertex_id_2,
                         stencilData[0] );
                  }
                  else
                  {
                     VertexDoFToEdgeDoF::generated::apply_3D_macroface_one_sided_vertexdof_to_edgedof_add(
                         &dstData[firstIdx[eo::X]],
                         &dstData[firstIdx[eo::XY]],
                         &dstData[firstIdx[eo::Y]],
                         &srcData[0],
                         &srcData[offset_gl_0],
                         static_cast< int32_t >( level ),
                         neighbor_cell_0_local_vertex_id_0,
                         neighbor_cell_0_local_vertex_id_1,
                         neighbor_cell_0_local_vertex_id_2,
                         stencilData[0] );
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

                     VertexDoFToEdgeDoF::generated::apply_3D_macroface_one_sided_vertexdof_to_edgedof_add(
                         &dstData[firstIdx[eo::X]],
                         &dstData[firstIdx[eo::XY]],
                         &dstData[firstIdx[eo::Y]],
                         &srcData[0],
                         &srcData[offset_gl_1],
                         static_cast< int32_t >( level ),
                         neighbor_cell_1_local_vertex_id_0,
                         neighbor_cell_1_local_vertex_id_1,
                         neighbor_cell_1_local_vertex_id_2,
                         stencilData[1] );
                  }
               }
               else
               {
                  VertexDoFToEdgeDoF::applyFace3D(
                      level, face, *storage_, faceStencil3DID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
               }
            }
            else
            {
               if ( hyteg::globalDefines::useGeneratedKernels )
               {
                  real_t* opr_data                      = face.getData( faceStencilID_ )->getPointer( level );
                  real_t* vertexToDiagonalEdgeStencil   = &opr_data[4];
                  real_t* vertexToHorizontalEdgeStencil = &opr_data[0];
                  real_t* vertexToVerticalEdgeStencil   = &opr_data[8];
                  real_t* src_data                      = face.getData( src.getFaceDataID() )->getPointer( level );
                  real_t* dst_data                      = face.getData( dst.getFaceDataID() )->getPointer( level );

                  typedef edgedof::EdgeDoFOrientation eo;
                  std::map< eo, uint_t >              firstIdx;
                  for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
                     firstIdx[e] = edgedof::macroface::index( level, 0, 0, e );

                  if ( updateType == hyteg::Replace )
                  {
                     VertexDoFToEdgeDoF::generated::apply_2D_macroface_vertexdof_to_edgedof_replace(
                         &dst_data[firstIdx[eo::X]],
                         &dst_data[firstIdx[eo::XY]],
                         &dst_data[firstIdx[eo::Y]],
                         src_data,
                         vertexToDiagonalEdgeStencil,
                         vertexToHorizontalEdgeStencil,
                         vertexToVerticalEdgeStencil,
                         static_cast< int32_t >( level ) );
                  }
                  else if ( updateType == hyteg::Add )
                  {
                     VertexDoFToEdgeDoF::generated::apply_2D_macroface_vertexdof_to_edgedof_add(
                         &dst_data[firstIdx[eo::X]],
                         &dst_data[firstIdx[eo::XY]],
                         &dst_data[firstIdx[eo::Y]],
                         src_data,
                         vertexToDiagonalEdgeStencil,
                         vertexToHorizontalEdgeStencil,
                         vertexToVerticalEdgeStencil,
                         static_cast< int32_t >( level ) );
                  }
               }
               else
               {
                  VertexDoFToEdgeDoF::applyFace(
                      level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
               }
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Face" );

   src.endCommunication< Face, Edge >( level );

   this->timingTree_->start( "Macro-Edge" );

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
            VertexDoFToEdgeDoF::applyEdge3D(
                level, edge, *getStorage(), edgeStencil3DID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
         }
         else
         {
            VertexDoFToEdgeDoF::applyEdge( level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
         }
      }
   }

   this->timingTree_->stop( "Macro-Edge" );

   this->stopTiming( "Apply" );
}

namespace VertexDoFToEdgeDoF {

uint_t macroEdgeVertexDoFToEdgeDoFStencilSize( const uint_t& level, const Primitive& primitive )
{
   WALBERLA_UNUSED( level );
   return 2 + primitive.getNumNeighborFaces();
}

uint_t macroFaceVertexDoFToEdgeDoFStencilSize( const uint_t& level, const Primitive& primitive )
{
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( primitive );
   return 4 + 4 + 4;
}

uint_t macroCellVertexDoFToEdgeDoFStencilSize( const uint_t& level, const Primitive& primitive )
{
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( primitive );
   return 7 * 27;
}
} // namespace VertexDoFToEdgeDoF

template class VertexDoFToEdgeDoFOperator< P2FenicsForm< hyteg::fenics::NoAssemble, hyteg::fenics::NoAssemble > >;
template class VertexDoFToEdgeDoFOperator<
    P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >;
template class VertexDoFToEdgeDoFOperator<
    P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >;

template class VertexDoFToEdgeDoFOperator<
    P2FenicsForm< p2_divt_cell_integral_0_otherwise, p2_tet_divt_tet_cell_integral_0_otherwise > >;
template class VertexDoFToEdgeDoFOperator<
    P2FenicsForm< p2_divt_cell_integral_1_otherwise, p2_tet_divt_tet_cell_integral_1_otherwise > >;
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble, p2_tet_divt_tet_cell_integral_2_otherwise > >;
template class VertexDoFToEdgeDoFOperator<
    P2FenicsForm< p2_div_cell_integral_0_otherwise, p2_tet_div_tet_cell_integral_0_otherwise > >;
template class VertexDoFToEdgeDoFOperator<
    P2FenicsForm< p2_div_cell_integral_1_otherwise, p2_tet_div_tet_cell_integral_1_otherwise > >;
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble, p2_tet_div_tet_cell_integral_2_otherwise > >;

template class VertexDoFToEdgeDoFOperator<
    P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_0_otherwise, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > >;
template class VertexDoFToEdgeDoFOperator<
    P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_1_otherwise, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > >;
template class VertexDoFToEdgeDoFOperator<
    P1ToP2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > >;

template class VertexDoFToEdgeDoFOperator<
    P2FenicsForm< p2_pspg_cell_integral_0_otherwise, p2_tet_pspg_tet_cell_integral_0_otherwise > >;

template class VertexDoFToEdgeDoFOperator< P2LinearCombinationForm >;
template class VertexDoFToEdgeDoFOperator< P2RowSumForm >;

// The following instantiations are required as building blocks in the P2ConstantEpsilon operator class
// clang-format off
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_0_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_0_otherwise > >;
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_1_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_1_otherwise > >;
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_2_otherwise > >;

template class VertexDoFToEdgeDoFOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_2_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_3_otherwise > >;
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_3_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_4_otherwise > >;
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_5_otherwise > >;

template class VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_6_otherwise > >;
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_7_otherwise > >;
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_8_otherwise > >;
// clang-format on

// The following instantiations are required as building blocks in the P2ConstantFullViscousOperator class
// clang-format off
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< p2_stokes_full_cell_integral_0_otherwise, p2_tet_stokes_full_tet_cell_integral_0_otherwise > >;
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< p2_stokes_full_cell_integral_1_otherwise, p2_tet_stokes_full_tet_cell_integral_1_otherwise > >;
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_2_otherwise > >;

template class VertexDoFToEdgeDoFOperator< P2FenicsForm< p2_stokes_full_cell_integral_2_otherwise, p2_tet_stokes_full_tet_cell_integral_3_otherwise > >;
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< p2_stokes_full_cell_integral_3_otherwise, p2_tet_stokes_full_tet_cell_integral_4_otherwise > >;
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_5_otherwise > >;

template class VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_6_otherwise > >;
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_7_otherwise > >;
template class VertexDoFToEdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_8_otherwise > >;
// clang-format on

} // namespace hyteg
