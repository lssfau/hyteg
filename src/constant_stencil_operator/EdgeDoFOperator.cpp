/*
 * Copyright (c) 2017-2022 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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
#include "EdgeDoFOperator.hpp"

#include "core/OpenMP.h"

#include "hyteg/HytegDefinitions.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroCell.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroEdge.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/edgedofspace/EdgeDoFPetsc.hpp"
#include "constant_stencil_operator/EdgeDoFGeneratedKernels/apply_2D_macroface_edgedof_to_edgedof_add.hpp"
#include "constant_stencil_operator/EdgeDoFGeneratedKernels/apply_2D_macroface_edgedof_to_edgedof_replace.hpp"
#include "constant_stencil_operator/EdgeDoFGeneratedKernels/apply_3D_macrocell_edgedof_to_edgedof_add.hpp"
#include "constant_stencil_operator/EdgeDoFGeneratedKernels/apply_3D_macrocell_edgedof_to_edgedof_replace.hpp"
#include "constant_stencil_operator/EdgeDoFGeneratedKernels/apply_3D_macroface_one_sided_edgedof_to_edgedof_add.hpp"
#include "constant_stencil_operator/EdgeDoFGeneratedKernels/apply_3D_macroface_one_sided_edgedof_to_edgedof_replace.hpp"
#include "hyteg/forms/P2LinearCombinationForm.hpp"
#include "hyteg/forms/P2RowSumForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/p2functionspace/variablestencil/P2VariableStencilCommon.hpp"



namespace hyteg {

using walberla::int_c;

template < class EdgeDoFForm >
EdgeDoFOperator< EdgeDoFForm >::EdgeDoFOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                 const uint_t                               minLevel,
                                                 const uint_t                               maxLevel )
: EdgeDoFOperator< EdgeDoFForm >( storage, minLevel, maxLevel, EdgeDoFForm() )
{}

template < class EdgeDoFForm >
EdgeDoFOperator< EdgeDoFForm >::EdgeDoFOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                 const uint_t                               minLevel,
                                                 const uint_t                               maxLevel,
                                                 const EdgeDoFForm&                         form )
: Operator( storage, minLevel, maxLevel )
, form_( form )
{
   auto edgeDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Edge > >(
       minLevel_, maxLevel_, macroEdgeEdgeDoFToEdgeDoFStencilSize );

   auto edge3DDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge > >( minLevel_,
                                                                                                                     maxLevel_ );

   auto faceDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Face > >(
       minLevel_, maxLevel_, macroFaceEdgeDoFToEdgeDoFStencilSize );

   auto face3DDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face > >( minLevel_,
                                                                                                                     maxLevel_ );

   auto cellDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell > >( minLevel_,
                                                                                                                     maxLevel_ );

   storage->addEdgeData( edgeStencilID_, edgeDataHandling, "VertexDoFToEdgeDoFOperatorEdgeStencil" );
   storage->addEdgeData( edgeStencil3DID_, edge3DDataHandling, "VertexDoFToEdgeDoFOperatorEdge3DStencil" );
   storage->addFaceData( faceStencilID_, faceDataHandling, "VertexDoFToEdgeDoFOperatorFaceStencil" );
   storage->addFaceData( faceStencil3DID_, face3DDataHandling, "VertexDoFToEdgeDoFOperatorFace3DStencil" );
   storage->addCellData( cellStencilID_, cellDataHandling, "VertexDoFToEdgeDoFOperatorCellStencil" );

   if ( this->getStorage()->hasGlobalCells() )
   {
      assembleEdgeToEdgeStencils< EdgeDoFForm >(
          storage, minLevel, maxLevel, edgeStencil3DID_, faceStencil3DID_, cellStencilID_, form_ );
   }
   else
   {
      assembleStencils();
   }
}

template < class EdgeDoFForm >
void EdgeDoFOperator< EdgeDoFForm >::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                               const EdgeDoFFunction< idx_t >&             src,
                                               const EdgeDoFFunction< idx_t >&             dst,
                                               size_t                                      level,
                                               DoFType                                     flag ) const
{
   auto storage = src.getStorage();

   for ( auto& it : this->getStorage()->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         if ( storage->hasGlobalCells() )
         {
            saveEdgeOperator3D(
                level, edge, *storage, this->getEdgeStencil3DID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat );
         }
         else
         {
            saveEdgeOperator( level, edge, this->getEdgeStencilID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat );
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
               saveFaceOperator3D(
                   level, face, *storage, this->getFaceStencil3DID(), src.getFaceDataID(), dst.getFaceDataID(), mat );
            }
            else
            {
               saveFaceOperator( level, face, this->getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat );
            }
         }
      }

      for ( auto& it : this->getStorage()->getCells() )
      {
         Cell& cell = *it.second;

         const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
         if ( testFlag( cellBC, flag ) )
         {
            saveCellOperator( level, cell, this->getCellStencilID(), src.getCellDataID(), dst.getCellDataID(), mat );
         }
      }
   }
}

template < class EdgeDoFForm >
void EdgeDoFOperator< EdgeDoFForm >::assembleStencils()
{
   using namespace P2Elements;

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
         while ( edgeIt->y() == 0 || edgeIt->x() == 0 ||
                 edgeIt->x() + edgeIt->y() == idx_t( hyteg::levelinfo::num_microedges_per_edge( level ) - 1 ) )
         {
            edgeIt++;
         }

         const Point3D horizontalMicroEdgePosition =
             faceBottomLeftCoords + ( walberla::real_c( edgeIt->x() * 2 + 1 ) * horizontalMicroEdgeOffset +
                                      walberla::real_c( edgeIt->y() * 2 ) * verticalMicroEdgeOffset );
         const Point3D verticalMicroEdgePosition =
             faceBottomLeftCoords + ( walberla::real_c( edgeIt->x() * 2 ) * horizontalMicroEdgeOffset +
                                      walberla::real_c( edgeIt->y() * 2 + 1 ) * verticalMicroEdgeOffset );
         const Point3D diagonalMicroEdgePosition = horizontalMicroEdgePosition + verticalMicroEdgeOffset;

         P2::variablestencil::assembleEdgeToEdgeStencil( form_,
                                                         { horizontalMicroEdgePosition + dirHO_W,
                                                           horizontalMicroEdgePosition + dirHO_E,
                                                           horizontalMicroEdgePosition + dirHO_NW },
                                                         { edgedof::stencilIndexFromHorizontalEdge( SD::EDGE_DI_N ),
                                                           edgedof::stencilIndexFromHorizontalEdge( SD::EDGE_VE_NW ),
                                                           edgedof::stencilIndexFromHorizontalEdge( SD::EDGE_HO_C ) },
                                                         vStencil );
         P2::variablestencil::assembleEdgeToEdgeStencil( form_,
                                                         { horizontalMicroEdgePosition + dirHO_W,
                                                           horizontalMicroEdgePosition + dirHO_E,
                                                           horizontalMicroEdgePosition + dirHO_SE },
                                                         { edgedof::stencilIndexFromHorizontalEdge( SD::EDGE_VE_SE ),
                                                           edgedof::stencilIndexFromHorizontalEdge( SD::EDGE_DI_S ),
                                                           edgedof::stencilIndexFromHorizontalEdge( SD::EDGE_HO_C ) },
                                                         vStencil );

         P2::variablestencil::assembleEdgeToEdgeStencil(
             form_,
             { verticalMicroEdgePosition + dirVE_N, verticalMicroEdgePosition + dirVE_S, verticalMicroEdgePosition + dirVE_NW },
             { edgedof::stencilIndexFromVerticalEdge( SD::EDGE_DI_W ),
               edgedof::stencilIndexFromVerticalEdge( SD::EDGE_HO_NW ),
               edgedof::stencilIndexFromVerticalEdge( SD::EDGE_VE_C ) },
             vStencil );
         P2::variablestencil::assembleEdgeToEdgeStencil(
             form_,
             { verticalMicroEdgePosition + dirVE_N, verticalMicroEdgePosition + dirVE_S, verticalMicroEdgePosition + dirVE_SE },
             { edgedof::stencilIndexFromVerticalEdge( SD::EDGE_HO_SE ),
               edgedof::stencilIndexFromVerticalEdge( SD::EDGE_DI_E ),
               edgedof::stencilIndexFromVerticalEdge( SD::EDGE_VE_C ) },
             vStencil );

         P2::variablestencil::assembleEdgeToEdgeStencil(
             form_,
             { diagonalMicroEdgePosition + dirDI_NW, diagonalMicroEdgePosition + dirDI_SE, diagonalMicroEdgePosition + dirDI_SW },
             { edgedof::stencilIndexFromDiagonalEdge( SD::EDGE_HO_S ),
               edgedof::stencilIndexFromDiagonalEdge( SD::EDGE_VE_W ),
               edgedof::stencilIndexFromDiagonalEdge( SD::EDGE_DI_C ) },
             vStencil );
         P2::variablestencil::assembleEdgeToEdgeStencil(
             form_,
             { diagonalMicroEdgePosition + dirDI_NW, diagonalMicroEdgePosition + dirDI_SE, diagonalMicroEdgePosition + dirDI_NE },
             { edgedof::stencilIndexFromDiagonalEdge( SD::EDGE_VE_E ),
               edgedof::stencilIndexFromDiagonalEdge( SD::EDGE_HO_N ),
               edgedof::stencilIndexFromDiagonalEdge( SD::EDGE_DI_C ) },
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

         real_t h = real_c( 1.0 ) / ( walberla::real_c( rowsize ) );

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

         horizontalMicroEdgePosition = leftCoords + real_c( 0.5 ) * dS_se;

         P2::variablestencil::assembleEdgeToEdgeStencil(
             form_,
             { horizontalMicroEdgePosition + dir_W, horizontalMicroEdgePosition + dir_E, horizontalMicroEdgePosition + dir_SE },
             { edgedof::stencilIndexFromHorizontalEdge( SD::EDGE_VE_SE ),
               edgedof::stencilIndexFromHorizontalEdge( SD::EDGE_DI_S ),
               edgedof::stencilIndexFromHorizontalEdge( SD::EDGE_HO_C ) },
             vStencil );

         if ( edge.getNumNeighborFaces() == 2 )
         {
            form_.setGeometryMap( faceN->getGeometryMap() );

            P2::variablestencil::assembleEdgeToEdgeStencil( form_,
                                                            { horizontalMicroEdgePosition + dir_W,
                                                              horizontalMicroEdgePosition + dir_E,
                                                              horizontalMicroEdgePosition + dir_NW },
                                                            { edgedof::stencilIndexFromHorizontalEdge( SD::EDGE_DI_N ),
                                                              edgedof::stencilIndexFromHorizontalEdge( SD::EDGE_VE_NW ),
                                                              edgedof::stencilIndexFromHorizontalEdge( SD::EDGE_HO_C ) },
                                                            vStencil );
         }
      }
   }
}

template < class EdgeDoFForm >
void EdgeDoFOperator< EdgeDoFForm >::apply( const EdgeDoFFunction< real_t >& src,
                                            const EdgeDoFFunction< real_t >& dst,
                                            uint_t                           level,
                                            DoFType                          flag,
                                            UpdateType                       updateType ) const
{
   this->startTiming( "Apply" );

   src.startCommunication< Edge, Face >( level );
   src.endCommunication< Edge, Face >( level );
   src.startCommunication< Face, Cell >( level );
   src.endCommunication< Face, Cell >( level );
   src.communicate< Cell, Face >( level );
   src.startCommunication< Face, Edge >( level );

   this->timingTree_->start( "Macro-Cell" );

   if ( level >= 1 )
   {
      std::vector< PrimitiveID > cellIDs = this->getStorage()->getCellIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
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
#ifdef HYTEG_USE_GENERATED_KERNELS
               typedef edgedof::EdgeDoFOrientation eo;
               auto                                dstData     = cell.getData( dst.getCellDataID() )->getPointer( level );
               auto                                srcData     = cell.getData( src.getCellDataID() )->getPointer( level );
               auto                                stencilData = cell.getData( cellStencilID_ )->getData( level );
               std::map< eo, uint_t >              firstIdx;
               for ( auto e : edgedof::allEdgeDoFOrientations )
                  firstIdx[e] = edgedof::macrocell::index( level, 0, 0, 0, e );

               if ( updateType == Replace )
               {
                  edgedof::macrocell::generated::apply_3D_macrocell_edgedof_to_edgedof_replace( &dstData[firstIdx[eo::X]],
                                                                                                &dstData[firstIdx[eo::XY]],
                                                                                                &dstData[firstIdx[eo::XYZ]],
                                                                                                &dstData[firstIdx[eo::XZ]],
                                                                                                &dstData[firstIdx[eo::Y]],
                                                                                                &dstData[firstIdx[eo::YZ]],
                                                                                                &dstData[firstIdx[eo::Z]],
                                                                                                &srcData[firstIdx[eo::X]],
                                                                                                &srcData[firstIdx[eo::XY]],
                                                                                                &srcData[firstIdx[eo::XYZ]],
                                                                                                &srcData[firstIdx[eo::XZ]],
                                                                                                &srcData[firstIdx[eo::Y]],
                                                                                                &srcData[firstIdx[eo::YZ]],
                                                                                                &srcData[firstIdx[eo::Z]],
                                                                                                stencilData,
                                                                                                static_cast< int32_t >( level ) );
               }
               else
               {
                  edgedof::macrocell::generated::apply_3D_macrocell_edgedof_to_edgedof_add( &dstData[firstIdx[eo::X]],
                                                                                            &dstData[firstIdx[eo::XY]],
                                                                                            &dstData[firstIdx[eo::XYZ]],
                                                                                            &dstData[firstIdx[eo::XZ]],
                                                                                            &dstData[firstIdx[eo::Y]],
                                                                                            &dstData[firstIdx[eo::YZ]],
                                                                                            &dstData[firstIdx[eo::Z]],
                                                                                            &srcData[firstIdx[eo::X]],
                                                                                            &srcData[firstIdx[eo::XY]],
                                                                                            &srcData[firstIdx[eo::XYZ]],
                                                                                            &srcData[firstIdx[eo::XZ]],
                                                                                            &srcData[firstIdx[eo::Y]],
                                                                                            &srcData[firstIdx[eo::YZ]],
                                                                                            &srcData[firstIdx[eo::Z]],
                                                                                            stencilData,
                                                                                            static_cast< int32_t >( level ) );
               }
#endif
            }
            else
            {
               edgedof::macrocell::apply( level, cell, cellStencilID_, src.getCellDataID(), dst.getCellDataID(), updateType );
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Cell" );

   this->timingTree_->start( "Macro-Face" );

   if ( level >= 1 )
   {
      std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
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
               if ( hyteg::globalDefines::useGeneratedKernels && face.getNumNeighborCells() == 2 )
               {
#ifdef HYTEG_USE_GENERATED_KERNELS
                  auto opr_data = face.getData( faceStencil3DID_ )->getData( level );
                  auto src_data = face.getData( src.getFaceDataID() )->getPointer( level );
                  auto dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );

                  const uint_t offset_x  = edgedof::macroface::index( level, 0, 0, edgedof::EdgeDoFOrientation::X );
                  const uint_t offset_xy = edgedof::macroface::index( level, 0, 0, edgedof::EdgeDoFOrientation::XY );
                  const uint_t offset_y  = edgedof::macroface::index( level, 0, 0, edgedof::EdgeDoFOrientation::Y );

                  std::map< uint_t, std::map< edgedof::EdgeDoFOrientation, uint_t > > offset_gl_orientation;
                  for ( uint_t gl = 0; gl < 2; gl++ )
                  {
                     for ( const auto& eo : edgedof::allEdgeDoFOrientations )
                     {
                        offset_gl_orientation[gl][eo] = edgedof::macroface::index( level, 0, 0, eo, gl );
                     }
                  }

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

                  if ( updateType == Replace )
                  {
                     edgedof::macroface::generated::apply_3D_macroface_one_sided_edgedof_to_edgedof_replace(
                         &dst_data[offset_x],
                         &dst_data[offset_xy],
                         &dst_data[offset_y],
                         &src_data[offset_x],
                         opr_data[0],
                         static_cast< int32_t >( level ),
                         neighbor_cell_0_local_vertex_id_0,
                         neighbor_cell_0_local_vertex_id_1,
                         neighbor_cell_0_local_vertex_id_2 );
                  }
                  else
                  {
                     edgedof::macroface::generated::apply_3D_macroface_one_sided_edgedof_to_edgedof_add(
                         &dst_data[offset_x],
                         &dst_data[offset_xy],
                         &dst_data[offset_y],
                         &src_data[offset_x],
                         &src_data[offset_xy],
                         &src_data[offset_y],
                         &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
                         &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
                         &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
                         &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
                         &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
                         &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
                         &src_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
                         opr_data[0],
                         static_cast< int32_t >( level ),
                         neighbor_cell_0_local_vertex_id_0,
                         neighbor_cell_0_local_vertex_id_1,
                         neighbor_cell_0_local_vertex_id_2 );
                  }

                  edgedof::macroface::generated::apply_3D_macroface_one_sided_edgedof_to_edgedof_add(
                      &dst_data[offset_x],
                      &dst_data[offset_xy],
                      &dst_data[offset_y],
                      &src_data[offset_x],
                      &src_data[offset_xy],
                      &src_data[offset_y],
                      &src_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::X]],
                      &src_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XY]],
                      &src_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XYZ]],
                      &src_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XZ]],
                      &src_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Y]],
                      &src_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::YZ]],
                      &src_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Z]],
                      opr_data[1],
                      static_cast< int32_t >( level ),
                      neighbor_cell_1_local_vertex_id_0,
                      neighbor_cell_1_local_vertex_id_1,
                      neighbor_cell_1_local_vertex_id_2 );
#endif
               }
               else
               {
                  edgedof::macroface::apply3D(
                      level, face, *storage_, faceStencil3DID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
               }
            }
            else
            {
               if ( hyteg::globalDefines::useGeneratedKernels )
               {
#ifdef HYTEG_USE_GENERATED_KERNELS
                  typedef edgedof::EdgeDoFOrientation eo;
                  real_t*                             opr_data = face.getData( faceStencilID_ )->getPointer( level );
                  real_t*                             src_data = face.getData( src.getFaceDataID() )->getPointer( level );
                  real_t*                             dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );
                  std::map< eo, uint_t >              firstIdx;
                  for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
                     firstIdx[e] = edgedof::macroface::index( level, 0, 0, e );

                  if ( updateType == hyteg::Replace )
                  {
                     edgedof::macroface::generated::apply_2D_macroface_edgedof_to_edgedof_replace(
                         &dst_data[firstIdx[eo::X]],
                         &dst_data[firstIdx[eo::XY]],
                         &dst_data[firstIdx[eo::Y]],
                         &src_data[firstIdx[eo::X]],
                         &src_data[firstIdx[eo::XY]],
                         &src_data[firstIdx[eo::Y]],
                         &opr_data[5],
                         &opr_data[0],
                         &opr_data[10],
                         static_cast< int32_t >( level ) );
                  }
                  else if ( updateType == hyteg::Add )
                  {
                     edgedof::macroface::generated::apply_2D_macroface_edgedof_to_edgedof_add( &dst_data[firstIdx[eo::X]],
                                                                                               &dst_data[firstIdx[eo::XY]],
                                                                                               &dst_data[firstIdx[eo::Y]],
                                                                                               &src_data[firstIdx[eo::X]],
                                                                                               &src_data[firstIdx[eo::XY]],
                                                                                               &src_data[firstIdx[eo::Y]],
                                                                                               &opr_data[5],
                                                                                               &opr_data[0],
                                                                                               &opr_data[10],
                                                                                               static_cast< int32_t >( level ) );
                  }
#endif
               }
               else
               {
                  edgedof::macroface::apply( level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
               }
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Face" );

   src.endCommunication< Face, Edge >( level );

   this->timingTree_->start( "Macro-Edge" );

   std::vector< PrimitiveID > edgeIDs = this->getStorage()->getEdgeIDs();
#ifdef HYTEG_BUILD_WITH_OPENMP
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
            edgedof::macroedge::apply3D(
                level, edge, *storage_, edgeStencil3DID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
         }
         else
         {
            edgedof::macroedge::apply( level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
         }
      }
   }

   this->timingTree_->stop( "Macro-Edge" );

   this->stopTiming( "Apply" );
}

template < class EdgeDoFForm >
const PrimitiveDataID< StencilMemory< real_t >, Edge >& EdgeDoFOperator< EdgeDoFForm >::getEdgeStencilID() const
{
   return edgeStencilID_;
}

template < class EdgeDoFForm >
const PrimitiveDataID< StencilMemory< real_t >, Face >& EdgeDoFOperator< EdgeDoFForm >::getFaceStencilID() const
{
   return faceStencilID_;
}

template < class EdgeDoFForm >
const PrimitiveDataID< LevelWiseMemory< edgedof::macroedge::StencilMap_T >, Edge >&
    EdgeDoFOperator< EdgeDoFForm >::getEdgeStencil3DID() const
{
   return edgeStencil3DID_;
}

template < class EdgeDoFForm >
const PrimitiveDataID< LevelWiseMemory< edgedof::macroface::StencilMap_T >, Face >&
    EdgeDoFOperator< EdgeDoFForm >::getFaceStencil3DID() const
{
   return faceStencil3DID_;
}

template < class EdgeDoFForm >
const PrimitiveDataID< LevelWiseMemory< edgedof::macrocell::StencilMap_T >, Cell >&
    EdgeDoFOperator< EdgeDoFForm >::getCellStencilID() const
{
   return cellStencilID_;
}

/// on edges only one stencil is required since only the horizontal edge DoFs belong to the edge
uint_t macroEdgeEdgeDoFToEdgeDoFStencilSize( const uint_t& level, const Primitive& primitive )
{
   WALBERLA_UNUSED( level );
   return 1 + 2 * primitive.getNumNeighborFaces();
}

/// on face three stencils are needed for horizontal, vertical and diagonal DoFs
uint_t macroFaceEdgeDoFToEdgeDoFStencilSize( const uint_t& level, const Primitive& primitive )
{
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( primitive );
   return 5 + 5 + 5;
}

uint_t macroCellEdgeDoFToEdgeDoFStencilSize( const uint_t& level, const Primitive& primitive )
{
   WALBERLA_UNUSED( level );
   WALBERLA_UNUSED( primitive );
   return 7 * 7 * 27;
}

template class EdgeDoFOperator< P2FenicsForm< hyteg::fenics::NoAssemble, hyteg::fenics::NoAssemble > >;
template class EdgeDoFOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >;
template class EdgeDoFOperator<
    P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >;

template class EdgeDoFOperator< P2FenicsForm< p2_divt_cell_integral_0_otherwise, p2_tet_divt_tet_cell_integral_0_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< p2_divt_cell_integral_1_otherwise, p2_tet_divt_tet_cell_integral_1_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< fenics::NoAssemble, p2_tet_divt_tet_cell_integral_2_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< p2_div_cell_integral_0_otherwise, p2_tet_div_tet_cell_integral_0_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< p2_div_cell_integral_1_otherwise, p2_tet_div_tet_cell_integral_1_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< fenics::NoAssemble, p2_tet_div_tet_cell_integral_2_otherwise > >;

template class EdgeDoFOperator< P2FenicsForm< p2_pspg_cell_integral_0_otherwise, p2_tet_pspg_tet_cell_integral_0_otherwise > >;

template class EdgeDoFOperator< P2LinearCombinationForm >;
template class EdgeDoFOperator< P2RowSumForm >;

// The following instantiations are required as building blocks in the P2ConstantEpsilon operator class
// clang-format off
template class EdgeDoFOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_0_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_0_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_1_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_1_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_2_otherwise > >;

template class EdgeDoFOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_2_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_3_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_3_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_4_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_5_otherwise > >;

template class EdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_6_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_7_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_8_otherwise > >;
// clang-format on

// The following instantiations are required as building blocks in the P2ConstantFullViscousOperator class
// clang-format off
template class EdgeDoFOperator< P2FenicsForm< p2_stokes_full_cell_integral_0_otherwise, p2_tet_stokes_full_tet_cell_integral_0_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< p2_stokes_full_cell_integral_1_otherwise, p2_tet_stokes_full_tet_cell_integral_1_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_2_otherwise > >;

template class EdgeDoFOperator< P2FenicsForm< p2_stokes_full_cell_integral_2_otherwise, p2_tet_stokes_full_tet_cell_integral_3_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< p2_stokes_full_cell_integral_3_otherwise, p2_tet_stokes_full_tet_cell_integral_4_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_5_otherwise > >;

template class EdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_6_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_7_otherwise > >;
template class EdgeDoFOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_8_otherwise > >;
// clang-format on
} // namespace hyteg
