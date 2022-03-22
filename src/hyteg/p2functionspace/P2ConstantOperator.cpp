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
#include "P2ConstantOperator.hpp"

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wfloat-conversion"
#pragma clang diagnostic ignored "-Wsign-conversion"
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wfloat-conversion"
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include "core/OpenMP.h"

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/forms/P2LinearCombinationForm.hpp"
#include "hyteg/forms/P2RowSumForm.hpp"
#include "hyteg/p2functionspace/P2Elements.hpp"
#include "hyteg/p2functionspace/P2MacroCell.hpp"
#include "hyteg/p2functionspace/P2MacroEdge.hpp"
#include "hyteg/p2functionspace/P2MacroFace.hpp"
#include "hyteg/p2functionspace/P2MacroVertex.hpp"
#include "hyteg/p2functionspace/generatedKernels/sor_2D_macroface_P2_update_edgedofs.hpp"
#include "hyteg/p2functionspace/generatedKernels/sor_2D_macroface_P2_update_vertexdofs.hpp"
#include "hyteg/p2functionspace/generatedKernels/sor_3D_macrocell_P2_update_edgedofs_by_type.hpp"
#include "hyteg/p2functionspace/generatedKernels/sor_3D_macrocell_P2_update_vertexdofs.hpp"
#include "hyteg/p2functionspace/generatedKernels/sor_3D_macrocell_P2_update_vertexdofs_backwards.hpp"
#include "hyteg/p2functionspace/generatedKernels/sor_3D_macroface_P2_update_edgedofs.hpp"
#include "hyteg/p2functionspace/generatedKernels/sor_3D_macroface_P2_update_edgedofs_backwards.hpp"
#include "hyteg/p2functionspace/generatedKernels/sor_3D_macroface_P2_update_edgedofs_one_sided.hpp"
#include "hyteg/p2functionspace/generatedKernels/sor_3D_macroface_P2_update_edgedofs_one_sided_backwards.hpp"
#include "hyteg/p2functionspace/generatedKernels/sor_3D_macroface_P2_update_vertexdofs.hpp"
#include "hyteg/p2functionspace/generatedKernels/sor_3D_macroface_P2_update_vertexdofs_backwards.hpp"
#include "hyteg/p2functionspace/generatedKernels/sor_3D_macroface_P2_update_vertexdofs_one_sided.hpp"
#include "hyteg/p2functionspace/generatedKernels/sor_3D_macroface_P2_update_vertexdofs_one_sided_backwards.hpp"

namespace hyteg {

using walberla::int_c;

template < class P2Form >
P2ConstantOperator< P2Form >::P2ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                  size_t                                     minLevel,
                                                  size_t                                     maxLevel )
: P2ConstantOperator< P2Form >( storage, minLevel, maxLevel, P2Form() )
{}

template < class P2Form >
P2ConstantOperator< P2Form >::P2ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                                  size_t                                     minLevel,
                                                  size_t                                     maxLevel,
                                                  const P2Form&                              form )
: Operator( storage, minLevel, maxLevel )
, vertexToVertex( storage, minLevel, maxLevel, form )
, edgeToVertex( storage, minLevel, maxLevel, form )
, vertexToEdge( storage, minLevel, maxLevel, form )
, edgeToEdge( storage, minLevel, maxLevel, form )
, form_( form )
{}

template < class P2Form >
void P2ConstantOperator< P2Form >::apply( const P2Function< real_t >& src,
                                          const P2Function< real_t >& dst,
                                          size_t                      level,
                                          DoFType                     flag,
                                          UpdateType                  updateType ) const
{
   WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

   if ( src.isDummy() )
   {
      return;
   }

   vertexToVertex.apply( src.getVertexDoFFunction(), dst.getVertexDoFFunction(), level, flag, updateType );
   edgeToVertex.apply( src.getEdgeDoFFunction(), dst.getVertexDoFFunction(), level, flag, Add );

   edgeToEdge.apply( src.getEdgeDoFFunction(), dst.getEdgeDoFFunction(), level, flag, updateType );
   vertexToEdge.apply( src.getVertexDoFFunction(), dst.getEdgeDoFFunction(), level, flag, Add );
}

template < class P2Form >
void P2ConstantOperator< P2Form >::smooth_gs( const P2Function< real_t >& dst,
                                              const P2Function< real_t >& rhs,
                                              const size_t                level,
                                              const DoFType               flag ) const
{
   smooth_sor( dst, rhs, 1.0, level, flag );
}

template < class P2Form >
void P2ConstantOperator< P2Form >::smooth_sor_macro_vertices( const P2Function< real_t >& dst,
                                                              const P2Function< real_t >& rhs,
                                                              const real_t&               relax,
                                                              const size_t                level,
                                                              const DoFType               flag,
                                                              const bool&                 backwards ) const
{
   WALBERLA_UNUSED( backwards );

   this->timingTree_->start( "Macro-Vertex" );

   for ( auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         if ( storage_->hasGlobalCells() )
         {
            P2::macrovertex::smoothSOR3D( level,
                                          *storage_,
                                          vertex,
                                          relax,
                                          vertexToVertex.getVertexStencilID(),
                                          edgeToVertex.getVertexStencil3DID(),
                                          dst.getVertexDoFFunction().getVertexDataID(),
                                          rhs.getVertexDoFFunction().getVertexDataID(),
                                          dst.getEdgeDoFFunction().getVertexDataID() );
         }
         else
         {
            P2::macrovertex::smoothSORVertexDoF( level,
                                                 vertex,
                                                 relax,
                                                 vertexToVertex.getVertexStencilID(),
                                                 dst.getVertexDoFFunction().getVertexDataID(),
                                                 edgeToVertex.getVertexStencilID(),
                                                 dst.getEdgeDoFFunction().getVertexDataID(),
                                                 rhs.getVertexDoFFunction().getVertexDataID() );
         }
      }
   }

   this->timingTree_->stop( "Macro-Vertex" );
}

template < class P2Form >
void P2ConstantOperator< P2Form >::smooth_sor_macro_edges( const P2Function< real_t >& dst,
                                                           const P2Function< real_t >& rhs,
                                                           const real_t&               relax,
                                                           const size_t                level,
                                                           const DoFType               flag,
                                                           const bool&                 backwards ) const
{
   this->timingTree_->start( "Macro-Edge" );

   std::vector< PrimitiveID::IDType > edgeIDs;
   for ( auto& it : storage_->getEdges() )
   {
      edgeIDs.push_back( it.first );
   }

#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
   {
      Edge& edge = *storage_->getEdge( PrimitiveID( edgeIDs[uint_c( i )] ) );

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         if ( storage_->hasGlobalCells() )
         {
            P2::macroedge::smoothSOR3D( level,
                                        *storage_,
                                        edge,
                                        relax,
                                        vertexToVertex.getEdgeStencilID(),
                                        vertexToVertex.getEdgeStencil3DID(),
                                        edgeToVertex.getEdgeStencil3DID(),
                                        vertexToEdge.getEdgeStencil3DID(),
                                        edgeToEdge.getEdgeStencil3DID(),
                                        dst.getVertexDoFFunction().getEdgeDataID(),
                                        rhs.getVertexDoFFunction().getEdgeDataID(),
                                        dst.getEdgeDoFFunction().getEdgeDataID(),
                                        rhs.getEdgeDoFFunction().getEdgeDataID(),
                                        backwards );
         }
         else
         {
            WALBERLA_CHECK( !backwards, "Backwards smoothing not implemented for P2 macroedge in 2D." )
            P2::macroedge::smoothSOR( level,
                                      edge,
                                      relax,
                                      vertexToVertex.getEdgeStencilID(),
                                      edgeToVertex.getEdgeStencilID(),
                                      dst.getVertexDoFFunction().getEdgeDataID(),
                                      vertexToEdge.getEdgeStencilID(),
                                      edgeToEdge.getEdgeStencilID(),
                                      dst.getEdgeDoFFunction().getEdgeDataID(),
                                      rhs.getVertexDoFFunction().getEdgeDataID(),
                                      rhs.getEdgeDoFFunction().getEdgeDataID() );
         }
      }
   }

   this->timingTree_->stop( "Macro-Edge" );
}

template < class P2Form >
void P2ConstantOperator< P2Form >::smooth_sor_macro_faces( const P2Function< real_t >& dst,
                                                           const P2Function< real_t >& rhs,
                                                           const real_t&               relax,
                                                           const size_t                level,
                                                           const DoFType               flag,
                                                           const bool&                 backwards ) const
{
   this->timingTree_->start( "Macro-Face" );

   std::vector< PrimitiveID::IDType > faceIDs;
   for ( auto& it : storage_->getFaces() )
   {
      faceIDs.push_back( it.first );
   }

#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( faceIDs.size() ); i++ )
   {
      Face& face = *storage_->getFace( PrimitiveID( faceIDs[uint_c( i )] ) );

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         if ( storage_->hasGlobalCells() )
         {
            if ( globalDefines::useGeneratedKernels )
            {
               using edgedof::EdgeDoFOrientation;
               using indexing::IndexIncrement;

               auto v2v_operator = face.getData( vertexToVertex.getFaceStencil3DID() )->getData( level );
               auto v2e_operator = face.getData( vertexToEdge.getFaceStencil3DID() )->getData( level );
               auto e2v_operator = face.getData( edgeToVertex.getFaceStencil3DID() )->getData( level );
               auto e2e_operator = face.getData( edgeToEdge.getFaceStencil3DID() )->getData( level );

               real_t* v_dst_data = face.getData( dst.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
               real_t* v_rhs_data = face.getData( rhs.getVertexDoFFunction().getFaceDataID() )->getPointer( level );

               real_t* e_dst_data = face.getData( dst.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
               real_t* e_rhs_data = face.getData( rhs.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

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

               if ( face.getNumNeighborCells() == 2 )
               {
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
                     if ( backwards )
                     {
                        P2::macroface::generated::sor_3D_macroface_P2_update_edgedofs_backwards(
                            &e_dst_data[offset_x],
                            &e_dst_data[offset_xy],
                            &e_dst_data[offset_y],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Z]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
                            &e_rhs_data[offset_x],
                            &e_rhs_data[offset_xy],
                            &e_rhs_data[offset_y],
                            v_dst_data,
                            &v_dst_data[vertex_offset_gl_1],
                            &v_dst_data[vertex_offset_gl_0],
                            e2e_operator[1],
                            e2e_operator[0],
                            static_cast< int32_t >( level ),
                            neighbor_cell_1_local_vertex_id_0,
                            neighbor_cell_1_local_vertex_id_1,
                            neighbor_cell_1_local_vertex_id_2,
                            neighbor_cell_0_local_vertex_id_0,
                            neighbor_cell_0_local_vertex_id_1,
                            neighbor_cell_0_local_vertex_id_2,
                            relax,
                            v2e_operator[1],
                            v2e_operator[0] );

                        P2::macroface::generated::sor_3D_macroface_P2_update_vertexdofs_backwards(
                            &e_dst_data[offset_x],
                            &e_dst_data[offset_xy],
                            &e_dst_data[offset_y],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Z]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
                            v_dst_data,
                            &v_dst_data[vertex_offset_gl_1],
                            &v_dst_data[vertex_offset_gl_0],
                            v_rhs_data,
                            e2v_operator[1],
                            e2v_operator[0],
                            static_cast< int32_t >( level ),
                            neighbor_cell_1_local_vertex_id_0,
                            neighbor_cell_1_local_vertex_id_1,
                            neighbor_cell_1_local_vertex_id_2,
                            neighbor_cell_0_local_vertex_id_0,
                            neighbor_cell_0_local_vertex_id_1,
                            neighbor_cell_0_local_vertex_id_2,
                            relax,
                            v2v_operator[1],
                            v2v_operator[0] );
                     }
                     else
                     {
                        P2::macroface::generated::sor_3D_macroface_P2_update_vertexdofs(
                            &e_dst_data[offset_x],
                            &e_dst_data[offset_xy],
                            &e_dst_data[offset_y],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Z]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
                            v_dst_data,
                            &v_dst_data[vertex_offset_gl_1],
                            &v_dst_data[vertex_offset_gl_0],
                            v_rhs_data,
                            e2v_operator[1],
                            e2v_operator[0],
                            static_cast< int32_t >( level ),
                            neighbor_cell_1_local_vertex_id_0,
                            neighbor_cell_1_local_vertex_id_1,
                            neighbor_cell_1_local_vertex_id_2,
                            neighbor_cell_0_local_vertex_id_0,
                            neighbor_cell_0_local_vertex_id_1,
                            neighbor_cell_0_local_vertex_id_2,
                            relax,
                            v2v_operator[1],
                            v2v_operator[0] );

                        P2::macroface::generated::sor_3D_macroface_P2_update_edgedofs(
                            &e_dst_data[offset_x],
                            &e_dst_data[offset_xy],
                            &e_dst_data[offset_y],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Z]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
                            &e_rhs_data[offset_x],
                            &e_rhs_data[offset_xy],
                            &e_rhs_data[offset_y],
                            v_dst_data,
                            &v_dst_data[vertex_offset_gl_1],
                            &v_dst_data[vertex_offset_gl_0],
                            e2e_operator[1],
                            e2e_operator[0],
                            static_cast< int32_t >( level ),
                            neighbor_cell_1_local_vertex_id_0,
                            neighbor_cell_1_local_vertex_id_1,
                            neighbor_cell_1_local_vertex_id_2,
                            neighbor_cell_0_local_vertex_id_0,
                            neighbor_cell_0_local_vertex_id_1,
                            neighbor_cell_0_local_vertex_id_2,
                            relax,
                            v2e_operator[1],
                            v2e_operator[0] );
                     }
                  }
                  else
                  {
                     if ( backwards )
                     {
                        P2::macroface::generated::sor_3D_macroface_P2_update_edgedofs_backwards(
                            &e_dst_data[offset_x],
                            &e_dst_data[offset_xy],
                            &e_dst_data[offset_y],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Z]],
                            &e_rhs_data[offset_x],
                            &e_rhs_data[offset_xy],
                            &e_rhs_data[offset_y],
                            v_dst_data,
                            &v_dst_data[vertex_offset_gl_0],
                            &v_dst_data[vertex_offset_gl_1],
                            e2e_operator[0],
                            e2e_operator[1],
                            static_cast< int32_t >( level ),
                            neighbor_cell_0_local_vertex_id_0,
                            neighbor_cell_0_local_vertex_id_1,
                            neighbor_cell_0_local_vertex_id_2,
                            neighbor_cell_1_local_vertex_id_0,
                            neighbor_cell_1_local_vertex_id_1,
                            neighbor_cell_1_local_vertex_id_2,
                            relax,
                            v2e_operator[0],
                            v2e_operator[1] );

                        P2::macroface::generated::sor_3D_macroface_P2_update_vertexdofs_backwards(
                            &e_dst_data[offset_x],
                            &e_dst_data[offset_xy],
                            &e_dst_data[offset_y],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Z]],
                            v_dst_data,
                            &v_dst_data[vertex_offset_gl_0],
                            &v_dst_data[vertex_offset_gl_1],
                            v_rhs_data,
                            e2v_operator[0],
                            e2v_operator[1],
                            static_cast< int32_t >( level ),
                            neighbor_cell_0_local_vertex_id_0,
                            neighbor_cell_0_local_vertex_id_1,
                            neighbor_cell_0_local_vertex_id_2,
                            neighbor_cell_1_local_vertex_id_0,
                            neighbor_cell_1_local_vertex_id_1,
                            neighbor_cell_1_local_vertex_id_2,
                            relax,
                            v2v_operator[0],
                            v2v_operator[1] );
                     }
                     else
                     {
                        P2::macroface::generated::sor_3D_macroface_P2_update_vertexdofs(
                            &e_dst_data[offset_x],
                            &e_dst_data[offset_xy],
                            &e_dst_data[offset_y],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Z]],
                            v_dst_data,
                            &v_dst_data[vertex_offset_gl_0],
                            &v_dst_data[vertex_offset_gl_1],
                            v_rhs_data,
                            e2v_operator[0],
                            e2v_operator[1],
                            static_cast< int32_t >( level ),
                            neighbor_cell_0_local_vertex_id_0,
                            neighbor_cell_0_local_vertex_id_1,
                            neighbor_cell_0_local_vertex_id_2,
                            neighbor_cell_1_local_vertex_id_0,
                            neighbor_cell_1_local_vertex_id_1,
                            neighbor_cell_1_local_vertex_id_2,
                            relax,
                            v2v_operator[0],
                            v2v_operator[1] );

                        P2::macroface::generated::sor_3D_macroface_P2_update_edgedofs(
                            &e_dst_data[offset_x],
                            &e_dst_data[offset_xy],
                            &e_dst_data[offset_y],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::X]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XY]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XYZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::XZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Y]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::YZ]],
                            &e_dst_data[offset_gl_orientation[1][edgedof::EdgeDoFOrientation::Z]],
                            &e_rhs_data[offset_x],
                            &e_rhs_data[offset_xy],
                            &e_rhs_data[offset_y],
                            v_dst_data,
                            &v_dst_data[vertex_offset_gl_0],
                            &v_dst_data[vertex_offset_gl_1],
                            e2e_operator[0],
                            e2e_operator[1],
                            static_cast< int32_t >( level ),
                            neighbor_cell_0_local_vertex_id_0,
                            neighbor_cell_0_local_vertex_id_1,
                            neighbor_cell_0_local_vertex_id_2,
                            neighbor_cell_1_local_vertex_id_0,
                            neighbor_cell_1_local_vertex_id_1,
                            neighbor_cell_1_local_vertex_id_2,
                            relax,
                            v2e_operator[0],
                            v2e_operator[1] );
                     }
                  }
               }
               else // only one neighbor face
               {
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

                  if ( backwards )
                  {
                     P2::macroface::generated::sor_3D_macroface_P2_update_edgedofs_one_sided_backwards(
                         &e_dst_data[offset_x],
                         &e_dst_data[offset_xy],
                         &e_dst_data[offset_y],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
                         &e_rhs_data[offset_x],
                         &e_rhs_data[offset_xy],
                         &e_rhs_data[offset_y],
                         v_dst_data,
                         &v_dst_data[vertex_offset_gl_0],
                         e2e_operator[0],
                         static_cast< int32_t >( level ),
                         neighbor_cell_0_local_vertex_id_0,
                         neighbor_cell_0_local_vertex_id_1,
                         neighbor_cell_0_local_vertex_id_2,
                         relax,
                         v2e_operator[0] );

                     P2::macroface::generated::sor_3D_macroface_P2_update_vertexdofs_one_sided_backwards(
                         &e_dst_data[offset_x],
                         &e_dst_data[offset_xy],
                         &e_dst_data[offset_y],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
                         v_dst_data,
                         &v_dst_data[vertex_offset_gl_0],
                         v_rhs_data,
                         e2v_operator[0],
                         static_cast< int32_t >( level ),
                         neighbor_cell_0_local_vertex_id_0,
                         neighbor_cell_0_local_vertex_id_1,
                         neighbor_cell_0_local_vertex_id_2,
                         relax,
                         v2v_operator[0] );
                  }
                  else
                  {
                     P2::macroface::generated::sor_3D_macroface_P2_update_vertexdofs_one_sided(
                         &e_dst_data[offset_x],
                         &e_dst_data[offset_xy],
                         &e_dst_data[offset_y],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
                         v_dst_data,
                         &v_dst_data[vertex_offset_gl_0],
                         v_rhs_data,
                         e2v_operator[0],
                         static_cast< int32_t >( level ),
                         neighbor_cell_0_local_vertex_id_0,
                         neighbor_cell_0_local_vertex_id_1,
                         neighbor_cell_0_local_vertex_id_2,
                         relax,
                         v2v_operator[0] );

                     P2::macroface::generated::sor_3D_macroface_P2_update_edgedofs_one_sided(
                         &e_dst_data[offset_x],
                         &e_dst_data[offset_xy],
                         &e_dst_data[offset_y],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::X]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XY]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XYZ]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::XZ]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Y]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::YZ]],
                         &e_dst_data[offset_gl_orientation[0][edgedof::EdgeDoFOrientation::Z]],
                         &e_rhs_data[offset_x],
                         &e_rhs_data[offset_xy],
                         &e_rhs_data[offset_y],
                         v_dst_data,
                         &v_dst_data[vertex_offset_gl_0],
                         e2e_operator[0],
                         static_cast< int32_t >( level ),
                         neighbor_cell_0_local_vertex_id_0,
                         neighbor_cell_0_local_vertex_id_1,
                         neighbor_cell_0_local_vertex_id_2,
                         relax,
                         v2e_operator[0] );
                  }
               }
            }
            else
            {
               WALBERLA_CHECK( !backwards );
               P2::macroface::smoothSOR3D( level,
                                           *storage_,
                                           face,
                                           relax,
                                           vertexToVertex.getFaceStencil3DID(),
                                           edgeToVertex.getFaceStencil3DID(),
                                           vertexToEdge.getFaceStencil3DID(),
                                           edgeToEdge.getFaceStencil3DID(),
                                           dst.getVertexDoFFunction().getFaceDataID(),
                                           rhs.getVertexDoFFunction().getFaceDataID(),
                                           dst.getEdgeDoFFunction().getFaceDataID(),
                                           rhs.getEdgeDoFFunction().getFaceDataID() );
            }
         }
         else
         {
            if ( globalDefines::useGeneratedKernels )
            {
               WALBERLA_CHECK( !backwards );

               real_t* v_dst_data = face.getData( dst.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
               real_t* v_rhs_data = face.getData( rhs.getVertexDoFFunction().getFaceDataID() )->getPointer( level );

               real_t* e_dst_data = face.getData( dst.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
               real_t* e_rhs_data = face.getData( rhs.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );

               real_t* v2v_opr_data = face.getData( vertexToVertex.getFaceStencilID() )->getPointer( level );
               real_t* v2e_opr_data = face.getData( vertexToEdge.getFaceStencilID() )->getPointer( level );
               real_t* e2v_opr_data = face.getData( edgeToVertex.getFaceStencilID() )->getPointer( level );
               real_t* e2e_opr_data = face.getData( edgeToEdge.getFaceStencilID() )->getPointer( level );

               typedef edgedof::EdgeDoFOrientation eo;
               std::map< eo, uint_t >              firstIdx;
               for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
                  firstIdx[e] = edgedof::macroface::index( level, 0, 0, e );

               P2::macroface::generated::sor_2D_macroface_P2_update_vertexdofs( &e_dst_data[firstIdx[eo::X]],
                                                                                &e_dst_data[firstIdx[eo::XY]],
                                                                                &e_dst_data[firstIdx[eo::Y]],
                                                                                e2v_opr_data,
                                                                                v_dst_data,
                                                                                v_rhs_data,
                                                                                v2v_opr_data,
                                                                                static_cast< int32_t >( level ),
                                                                                relax );
               P2::macroface::generated::sor_2D_macroface_P2_update_edgedofs( &e_dst_data[firstIdx[eo::X]],
                                                                              &e_dst_data[firstIdx[eo::XY]],
                                                                              &e_dst_data[firstIdx[eo::Y]],
                                                                              &e_rhs_data[firstIdx[eo::X]],
                                                                              &e_rhs_data[firstIdx[eo::XY]],
                                                                              &e_rhs_data[firstIdx[eo::Y]],
                                                                              &e2e_opr_data[0],
                                                                              &e2e_opr_data[5],
                                                                              &e2e_opr_data[10],
                                                                              v_dst_data,
                                                                              &v2e_opr_data[0],
                                                                              &v2e_opr_data[4],
                                                                              &v2e_opr_data[8],
                                                                              static_cast< int32_t >( level ),
                                                                              relax );
            }
            else
            {
               WALBERLA_CHECK( !backwards );
               P2::macroface::smoothSOR( level,
                                         face,
                                         relax,
                                         vertexToVertex.getFaceStencilID(),
                                         edgeToVertex.getFaceStencilID(),
                                         dst.getVertexDoFFunction().getFaceDataID(),
                                         vertexToEdge.getFaceStencilID(),
                                         edgeToEdge.getFaceStencilID(),
                                         dst.getEdgeDoFFunction().getFaceDataID(),
                                         rhs.getVertexDoFFunction().getFaceDataID(),
                                         rhs.getEdgeDoFFunction().getFaceDataID() );
            }
         }
      }
   }

   this->timingTree_->stop( "Macro-Face" );
}

template < class P2Form >
void P2ConstantOperator< P2Form >::smooth_sor_macro_cells( const P2Function< real_t >& dst,
                                                           const P2Function< real_t >& rhs,
                                                           const real_t&               relax,
                                                           const size_t                level,
                                                           const DoFType               flag,
                                                           const bool&                 backwards ) const
{
   this->timingTree_->start( "Macro-Cell" );

   std::vector< PrimitiveID::IDType > cellIDs;
   for ( auto& it : storage_->getCells() )
   {
      cellIDs.push_back( it.first );
   }

#ifdef WALBERLA_BUILD_WITH_OPENMP
#pragma omp parallel for default( shared )
#endif
   for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
   {
      Cell& cell = *storage_->getCell( PrimitiveID( cellIDs[uint_c( i )] ) );

      const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
      if ( testFlag( cellBC, flag ) )
      {
         if ( globalDefines::useGeneratedKernels )
         {
            real_t* v_dst_data = cell.getData( dst.getVertexDoFFunction().getCellDataID() )->getPointer( level );
            real_t* v_rhs_data = cell.getData( rhs.getVertexDoFFunction().getCellDataID() )->getPointer( level );

            real_t* e_dst_data = cell.getData( dst.getEdgeDoFFunction().getCellDataID() )->getPointer( level );
            real_t* e_rhs_data = cell.getData( rhs.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

            auto v2v_opr_data = cell.getData( vertexToVertex.getCellStencilID() )->getData( level );
            auto v2e_opr_data = cell.getData( vertexToEdge.getCellStencilID() )->getData( level );
            auto e2v_opr_data = cell.getData( edgeToVertex.getCellStencilID() )->getData( level );
            auto e2e_opr_data = cell.getData( edgeToEdge.getCellStencilID() )->getData( level );

            typedef edgedof::EdgeDoFOrientation eo;
            std::map< eo, uint_t >              firstIdx;
            for ( auto e : edgedof::allEdgeDoFOrientations )
               firstIdx[e] = edgedof::macrocell::index( level, 0, 0, 0, e );

            if ( backwards )
            {
               WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start( "Updating EdgeDoFs" ); }

               // Splitting the SOR into multiple sweeps: one per edge type.
               // This has severe performance advantages.
               // Due to the memory layout of the edge DoFs and the coincidence that there are
               // never two (or more) edges per type in one element, this split is a natural
               // coloring of the edge DoFs.
               // Therefore each of the split-kernels can be vectorized!
               // However it is not clear how the smoothing property suffers from this splitting.
               // In first tests it is marginally worse - but the performance gain is huge.

               P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_XYZ( &e_dst_data[firstIdx[eo::X]],
                                                                                          &e_dst_data[firstIdx[eo::XY]],
                                                                                          &e_dst_data[firstIdx[eo::XYZ]],
                                                                                          &e_dst_data[firstIdx[eo::XZ]],
                                                                                          &e_dst_data[firstIdx[eo::Y]],
                                                                                          &e_dst_data[firstIdx[eo::YZ]],
                                                                                          &e_dst_data[firstIdx[eo::Z]],
                                                                                          &e_rhs_data[firstIdx[eo::XYZ]],
                                                                                          v_dst_data,
                                                                                          e2e_opr_data,
                                                                                          static_cast< int32_t >( level ),
                                                                                          relax,
                                                                                          v2e_opr_data );

               P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_YZ( &e_dst_data[firstIdx[eo::X]],
                                                                                         &e_dst_data[firstIdx[eo::XY]],
                                                                                         &e_dst_data[firstIdx[eo::XYZ]],
                                                                                         &e_dst_data[firstIdx[eo::XZ]],
                                                                                         &e_dst_data[firstIdx[eo::Y]],
                                                                                         &e_dst_data[firstIdx[eo::YZ]],
                                                                                         &e_dst_data[firstIdx[eo::Z]],
                                                                                         &e_rhs_data[firstIdx[eo::YZ]],
                                                                                         v_dst_data,
                                                                                         e2e_opr_data,
                                                                                         static_cast< int32_t >( level ),
                                                                                         relax,
                                                                                         v2e_opr_data );

               P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_XZ( &e_dst_data[firstIdx[eo::X]],
                                                                                         &e_dst_data[firstIdx[eo::XY]],
                                                                                         &e_dst_data[firstIdx[eo::XYZ]],
                                                                                         &e_dst_data[firstIdx[eo::XZ]],
                                                                                         &e_dst_data[firstIdx[eo::Y]],
                                                                                         &e_dst_data[firstIdx[eo::YZ]],
                                                                                         &e_dst_data[firstIdx[eo::Z]],
                                                                                         &e_rhs_data[firstIdx[eo::XZ]],
                                                                                         v_dst_data,
                                                                                         e2e_opr_data,
                                                                                         static_cast< int32_t >( level ),
                                                                                         relax,
                                                                                         v2e_opr_data );

               P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_XY( &e_dst_data[firstIdx[eo::X]],
                                                                                         &e_dst_data[firstIdx[eo::XY]],
                                                                                         &e_dst_data[firstIdx[eo::XYZ]],
                                                                                         &e_dst_data[firstIdx[eo::XZ]],
                                                                                         &e_dst_data[firstIdx[eo::Y]],
                                                                                         &e_dst_data[firstIdx[eo::YZ]],
                                                                                         &e_dst_data[firstIdx[eo::Z]],
                                                                                         &e_rhs_data[firstIdx[eo::XY]],
                                                                                         v_dst_data,
                                                                                         e2e_opr_data,
                                                                                         static_cast< int32_t >( level ),
                                                                                         relax,
                                                                                         v2e_opr_data );

               P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_Z( &e_dst_data[firstIdx[eo::X]],
                                                                                        &e_dst_data[firstIdx[eo::XY]],
                                                                                        &e_dst_data[firstIdx[eo::XYZ]],
                                                                                        &e_dst_data[firstIdx[eo::XZ]],
                                                                                        &e_dst_data[firstIdx[eo::Y]],
                                                                                        &e_dst_data[firstIdx[eo::YZ]],
                                                                                        &e_dst_data[firstIdx[eo::Z]],
                                                                                        &e_rhs_data[firstIdx[eo::Z]],
                                                                                        v_dst_data,
                                                                                        e2e_opr_data,
                                                                                        static_cast< int32_t >( level ),
                                                                                        relax,
                                                                                        v2e_opr_data );

               P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_Y( &e_dst_data[firstIdx[eo::X]],
                                                                                        &e_dst_data[firstIdx[eo::XY]],
                                                                                        &e_dst_data[firstIdx[eo::XYZ]],
                                                                                        &e_dst_data[firstIdx[eo::XZ]],
                                                                                        &e_dst_data[firstIdx[eo::Y]],
                                                                                        &e_dst_data[firstIdx[eo::YZ]],
                                                                                        &e_dst_data[firstIdx[eo::Z]],
                                                                                        &e_rhs_data[firstIdx[eo::Y]],
                                                                                        v_dst_data,
                                                                                        e2e_opr_data,
                                                                                        static_cast< int32_t >( level ),
                                                                                        relax,
                                                                                        v2e_opr_data );

               P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_X( &e_dst_data[firstIdx[eo::X]],
                                                                                        &e_dst_data[firstIdx[eo::XY]],
                                                                                        &e_dst_data[firstIdx[eo::XYZ]],
                                                                                        &e_dst_data[firstIdx[eo::XZ]],
                                                                                        &e_dst_data[firstIdx[eo::Y]],
                                                                                        &e_dst_data[firstIdx[eo::YZ]],
                                                                                        &e_dst_data[firstIdx[eo::Z]],
                                                                                        &e_rhs_data[firstIdx[eo::X]],
                                                                                        v_dst_data,
                                                                                        e2e_opr_data,
                                                                                        static_cast< int32_t >( level ),
                                                                                        relax,
                                                                                        v2e_opr_data );

               WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->stop( "Updating EdgeDoFs" ); }

               WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start( "Updating VertexDoFs" ); }

               P2::macrocell::generated::sor_3D_macrocell_P2_update_vertexdofs_backwards( &e_dst_data[firstIdx[eo::X]],
                                                                                          &e_dst_data[firstIdx[eo::XY]],
                                                                                          &e_dst_data[firstIdx[eo::XYZ]],
                                                                                          &e_dst_data[firstIdx[eo::XZ]],
                                                                                          &e_dst_data[firstIdx[eo::Y]],
                                                                                          &e_dst_data[firstIdx[eo::YZ]],
                                                                                          &e_dst_data[firstIdx[eo::Z]],
                                                                                          v_dst_data,
                                                                                          v_rhs_data,
                                                                                          e2v_opr_data,
                                                                                          static_cast< int32_t >( level ),
                                                                                          relax,
                                                                                          v2v_opr_data );

               WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->stop( "Updating VertexDoFs" ); }
            }
            else
            {
               WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start( "Updating VertexDoFs" ); }

               P2::macrocell::generated::sor_3D_macrocell_P2_update_vertexdofs( &e_dst_data[firstIdx[eo::X]],
                                                                                &e_dst_data[firstIdx[eo::XY]],
                                                                                &e_dst_data[firstIdx[eo::XYZ]],
                                                                                &e_dst_data[firstIdx[eo::XZ]],
                                                                                &e_dst_data[firstIdx[eo::Y]],
                                                                                &e_dst_data[firstIdx[eo::YZ]],
                                                                                &e_dst_data[firstIdx[eo::Z]],
                                                                                v_dst_data,
                                                                                v_rhs_data,
                                                                                e2v_opr_data,
                                                                                static_cast< int32_t >( level ),
                                                                                relax,
                                                                                v2v_opr_data );

               WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->stop( "Updating VertexDoFs" ); }

               WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start( "Updating EdgeDoFs" ); }

               // Splitting the SOR into multiple sweeps: one per edge type.
               // This has severe performance advantages.
               // Due to the memory layout of the edge DoFs and the coincidence that there are
               // never two (or more) edges per type in one element, this split is a natural
               // coloring of the edge DoFs.
               // Therefore each of the split-kernels can be vectorized!
               // However it is not clear how the smoothing property suffers from this splitting.
               // In first tests it is marginally worse - but the performance gain is huge.

               P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_X( &e_dst_data[firstIdx[eo::X]],
                                                                                        &e_dst_data[firstIdx[eo::XY]],
                                                                                        &e_dst_data[firstIdx[eo::XYZ]],
                                                                                        &e_dst_data[firstIdx[eo::XZ]],
                                                                                        &e_dst_data[firstIdx[eo::Y]],
                                                                                        &e_dst_data[firstIdx[eo::YZ]],
                                                                                        &e_dst_data[firstIdx[eo::Z]],
                                                                                        &e_rhs_data[firstIdx[eo::X]],
                                                                                        v_dst_data,
                                                                                        e2e_opr_data,
                                                                                        static_cast< int32_t >( level ),
                                                                                        relax,
                                                                                        v2e_opr_data );

               P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_Y( &e_dst_data[firstIdx[eo::X]],
                                                                                        &e_dst_data[firstIdx[eo::XY]],
                                                                                        &e_dst_data[firstIdx[eo::XYZ]],
                                                                                        &e_dst_data[firstIdx[eo::XZ]],
                                                                                        &e_dst_data[firstIdx[eo::Y]],
                                                                                        &e_dst_data[firstIdx[eo::YZ]],
                                                                                        &e_dst_data[firstIdx[eo::Z]],
                                                                                        &e_rhs_data[firstIdx[eo::Y]],
                                                                                        v_dst_data,
                                                                                        e2e_opr_data,
                                                                                        static_cast< int32_t >( level ),
                                                                                        relax,
                                                                                        v2e_opr_data );

               P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_Z( &e_dst_data[firstIdx[eo::X]],
                                                                                        &e_dst_data[firstIdx[eo::XY]],
                                                                                        &e_dst_data[firstIdx[eo::XYZ]],
                                                                                        &e_dst_data[firstIdx[eo::XZ]],
                                                                                        &e_dst_data[firstIdx[eo::Y]],
                                                                                        &e_dst_data[firstIdx[eo::YZ]],
                                                                                        &e_dst_data[firstIdx[eo::Z]],
                                                                                        &e_rhs_data[firstIdx[eo::Z]],
                                                                                        v_dst_data,
                                                                                        e2e_opr_data,
                                                                                        static_cast< int32_t >( level ),
                                                                                        relax,
                                                                                        v2e_opr_data );

               P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_XY( &e_dst_data[firstIdx[eo::X]],
                                                                                         &e_dst_data[firstIdx[eo::XY]],
                                                                                         &e_dst_data[firstIdx[eo::XYZ]],
                                                                                         &e_dst_data[firstIdx[eo::XZ]],
                                                                                         &e_dst_data[firstIdx[eo::Y]],
                                                                                         &e_dst_data[firstIdx[eo::YZ]],
                                                                                         &e_dst_data[firstIdx[eo::Z]],
                                                                                         &e_rhs_data[firstIdx[eo::XY]],
                                                                                         v_dst_data,
                                                                                         e2e_opr_data,
                                                                                         static_cast< int32_t >( level ),
                                                                                         relax,
                                                                                         v2e_opr_data );

               P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_XZ( &e_dst_data[firstIdx[eo::X]],
                                                                                         &e_dst_data[firstIdx[eo::XY]],
                                                                                         &e_dst_data[firstIdx[eo::XYZ]],
                                                                                         &e_dst_data[firstIdx[eo::XZ]],
                                                                                         &e_dst_data[firstIdx[eo::Y]],
                                                                                         &e_dst_data[firstIdx[eo::YZ]],
                                                                                         &e_dst_data[firstIdx[eo::Z]],
                                                                                         &e_rhs_data[firstIdx[eo::XZ]],
                                                                                         v_dst_data,
                                                                                         e2e_opr_data,
                                                                                         static_cast< int32_t >( level ),
                                                                                         relax,
                                                                                         v2e_opr_data );

               P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_YZ( &e_dst_data[firstIdx[eo::X]],
                                                                                         &e_dst_data[firstIdx[eo::XY]],
                                                                                         &e_dst_data[firstIdx[eo::XYZ]],
                                                                                         &e_dst_data[firstIdx[eo::XZ]],
                                                                                         &e_dst_data[firstIdx[eo::Y]],
                                                                                         &e_dst_data[firstIdx[eo::YZ]],
                                                                                         &e_dst_data[firstIdx[eo::Z]],
                                                                                         &e_rhs_data[firstIdx[eo::YZ]],
                                                                                         v_dst_data,
                                                                                         e2e_opr_data,
                                                                                         static_cast< int32_t >( level ),
                                                                                         relax,
                                                                                         v2e_opr_data );

               P2::macrocell::generated::sor_3D_macrocell_P2_update_edgedofs_by_type_XYZ( &e_dst_data[firstIdx[eo::X]],
                                                                                          &e_dst_data[firstIdx[eo::XY]],
                                                                                          &e_dst_data[firstIdx[eo::XYZ]],
                                                                                          &e_dst_data[firstIdx[eo::XZ]],
                                                                                          &e_dst_data[firstIdx[eo::Y]],
                                                                                          &e_dst_data[firstIdx[eo::YZ]],
                                                                                          &e_dst_data[firstIdx[eo::Z]],
                                                                                          &e_rhs_data[firstIdx[eo::XYZ]],
                                                                                          v_dst_data,
                                                                                          e2e_opr_data,
                                                                                          static_cast< int32_t >( level ),
                                                                                          relax,
                                                                                          v2e_opr_data );

               WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->stop( "Updating EdgeDoFs" ); }
            }
         }
         else
         {
            WALBERLA_CHECK( !backwards );
            P2::macrocell::smoothSOR( level,
                                      cell,
                                      relax,
                                      vertexToVertex.getCellStencilID(),
                                      edgeToVertex.getCellStencilID(),
                                      vertexToEdge.getCellStencilID(),
                                      edgeToEdge.getCellStencilID(),
                                      dst.getVertexDoFFunction().getCellDataID(),
                                      rhs.getVertexDoFFunction().getCellDataID(),
                                      dst.getEdgeDoFFunction().getCellDataID(),
                                      rhs.getEdgeDoFFunction().getCellDataID() );
         }
      }
   }

   this->timingTree_->stop( "Macro-Cell" );
}

template < class P2Form >
void P2ConstantOperator< P2Form >::smooth_sor( const P2Function< real_t >& dst,
                                               const P2Function< real_t >& rhs,
                                               real_t                      relax,
                                               const size_t                level,
                                               const DoFType               flag,
                                               const bool&                 backwards ) const
{
   if ( backwards )
   {
      WALBERLA_CHECK( globalDefines::useGeneratedKernels, "Backward SOR only implemented in generated kernels." )
      WALBERLA_CHECK( storage_->hasGlobalCells(), "Backward SOR currently only implemented for 3D for P2." )
      this->startTiming( "SOR backwards" );
   }
   else
   {
      this->startTiming( "SOR" );
   }

   communication::syncP2FunctionBetweenPrimitives( dst, level );

   if ( backwards )
   {
      smooth_sor_macro_cells( dst, rhs, relax, level, flag, backwards );

      dst.getVertexDoFFunction().communicate< Cell, Face >( level );
      dst.getEdgeDoFFunction().communicate< Cell, Face >( level );

      smooth_sor_macro_faces( dst, rhs, relax, level, flag, backwards );

      dst.getVertexDoFFunction().communicate< Face, Edge >( level );
      dst.getEdgeDoFFunction().communicate< Face, Edge >( level );

      smooth_sor_macro_edges( dst, rhs, relax, level, flag, backwards );

      dst.getVertexDoFFunction().communicate< Edge, Vertex >( level );
      dst.getEdgeDoFFunction().communicate< Edge, Vertex >( level );

      smooth_sor_macro_vertices( dst, rhs, relax, level, flag, backwards );
   }
   else
   {
      smooth_sor_macro_vertices( dst, rhs, relax, level, flag, backwards );

      dst.getVertexDoFFunction().communicate< Vertex, Edge >( level );
      dst.getEdgeDoFFunction().communicate< Vertex, Edge >( level );

      smooth_sor_macro_edges( dst, rhs, relax, level, flag, backwards );

      dst.getVertexDoFFunction().communicate< Edge, Face >( level );
      dst.getEdgeDoFFunction().communicate< Edge, Face >( level );

      smooth_sor_macro_faces( dst, rhs, relax, level, flag, backwards );

      dst.getVertexDoFFunction().communicate< Face, Cell >( level );
      dst.getEdgeDoFFunction().communicate< Face, Cell >( level );

      smooth_sor_macro_cells( dst, rhs, relax, level, flag, backwards );
   }

   if ( backwards )
      this->stopTiming( "SOR backwards" );
   else
      this->stopTiming( "SOR" );
}

template < class P2Form >
void P2ConstantOperator< P2Form >::smooth_jac( const P2Function< real_t >& dst,
                                               const P2Function< real_t >& rhs,
                                               const P2Function< real_t >& src,
                                               real_t                      relax,
                                               size_t                      level,
                                               DoFType                     flag ) const
{
   ///TODO: remove unneccessary communication here
   // src.getVertexDoFFunction().communicate< Face, Edge >( level );
   // src.getVertexDoFFunction().communicate< Edge, Vertex >( level );
   // src.getVertexDoFFunction().communicate< Vertex, Edge >( level );
   // src.getVertexDoFFunction().communicate< Edge, Face >( level );
   // src.getEdgeDoFFunction().communicate< Face, Edge >( level );
   // src.getEdgeDoFFunction().communicate< Edge, Vertex >( level );
   // src.getEdgeDoFFunction().communicate< Vertex, Edge >( level );
   // src.getEdgeDoFFunction().communicate< Edge, Face >( level );
   communication::syncP2FunctionBetweenPrimitives( src, level );

   if ( storage_->hasGlobalCells() )
   {
      throw std::runtime_error( "P2ConstantOperator::smooth_jac() not implemented for 3D, yet!" );
   }
   else
   {
      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            P2::macroface::smoothJacobiVertexDoF( level,
                                                  face,
                                                  vertexToVertex.getFaceStencilID(),
                                                  src.getVertexDoFFunction().getFaceDataID(),
                                                  dst.getVertexDoFFunction().getFaceDataID(),
                                                  edgeToVertex.getFaceStencilID(),
                                                  src.getEdgeDoFFunction().getFaceDataID(),
                                                  rhs.getVertexDoFFunction().getFaceDataID(),
                                                  relax );
         }
      }
      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            P2::macroface::smoothJacobiEdgeDoF( level,
                                                face,
                                                vertexToEdge.getFaceStencilID(),
                                                src.getVertexDoFFunction().getFaceDataID(),
                                                edgeToEdge.getFaceStencilID(),
                                                src.getEdgeDoFFunction().getFaceDataID(),
                                                dst.getEdgeDoFFunction().getFaceDataID(),
                                                rhs.getEdgeDoFFunction().getFaceDataID(),
                                                relax );
         }
      }

      for ( auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            P2::macroedge::smoothJacobi( level,
                                         edge,
                                         relax,
                                         vertexToVertex.getEdgeStencilID(),
                                         edgeToVertex.getEdgeStencilID(),
                                         src.getVertexDoFFunction().getEdgeDataID(),
                                         dst.getVertexDoFFunction().getEdgeDataID(),
                                         vertexToEdge.getEdgeStencilID(),
                                         edgeToEdge.getEdgeStencilID(),
                                         src.getEdgeDoFFunction().getEdgeDataID(),
                                         dst.getEdgeDoFFunction().getEdgeDataID(),
                                         rhs.getVertexDoFFunction().getEdgeDataID(),
                                         rhs.getEdgeDoFFunction().getEdgeDataID() );
         }
      }

      for ( auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if ( testFlag( vertexBC, flag ) )
         {
            P2::macrovertex::smoothJacobiVertexDoF( level,
                                                    vertex,
                                                    relax,
                                                    vertexToVertex.getVertexStencilID(),
                                                    src.getVertexDoFFunction().getVertexDataID(),
                                                    dst.getVertexDoFFunction().getVertexDataID(),
                                                    edgeToVertex.getVertexStencilID(),
                                                    src.getEdgeDoFFunction().getVertexDataID(),
                                                    rhs.getVertexDoFFunction().getVertexDataID() );
         }
      }
   }
}

template class P2ConstantOperator<
    P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >;

template class P2ConstantOperator< P2FenicsForm< p2_divt_cell_integral_0_otherwise, p2_tet_divt_tet_cell_integral_0_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< p2_divt_cell_integral_1_otherwise, p2_tet_divt_tet_cell_integral_1_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< fenics::NoAssemble, p2_tet_divt_tet_cell_integral_2_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< p2_div_cell_integral_0_otherwise, p2_tet_div_tet_cell_integral_0_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< p2_div_cell_integral_1_otherwise, p2_tet_div_tet_cell_integral_1_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< fenics::NoAssemble, p2_tet_div_tet_cell_integral_2_otherwise > >;

template class P2ConstantOperator< P2FenicsForm< p2_pspg_cell_integral_0_otherwise, p2_tet_pspg_tet_cell_integral_0_otherwise > >;

template class P2ConstantOperator< P2LinearCombinationForm >;
template class P2ConstantOperator< P2RowSumForm >;

// The following instantiations are required as building blocks in the P2ConstantEpsilonOperator class
// clang-format off
template class P2ConstantOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_0_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_0_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_1_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_1_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_2_otherwise > >;

template class P2ConstantOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_2_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_3_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< p2_stokes_epsilon_cell_integral_3_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_4_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_5_otherwise > >;

template class P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_6_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_7_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_8_otherwise > >;
// clang-format on

// The following instantiations are required as building blocks in the P2ConstantFullViscousOperator class
// clang-format off
template class P2ConstantOperator< P2FenicsForm< p2_stokes_full_cell_integral_0_otherwise, p2_tet_stokes_full_tet_cell_integral_0_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< p2_stokes_full_cell_integral_1_otherwise, p2_tet_stokes_full_tet_cell_integral_1_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_2_otherwise > >;

template class P2ConstantOperator< P2FenicsForm< p2_stokes_full_cell_integral_2_otherwise, p2_tet_stokes_full_tet_cell_integral_3_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< p2_stokes_full_cell_integral_3_otherwise, p2_tet_stokes_full_tet_cell_integral_4_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_5_otherwise > >;

template class P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_6_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_7_otherwise > >;
template class P2ConstantOperator< P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_8_otherwise > >;
// clang-format on

} // namespace hyteg
