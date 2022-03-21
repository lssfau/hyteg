/*
 * Copyright (c) 2017-2020 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl, Benjamin Mann.
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
#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include "hyteg/p1functionspace/P1ConstantOperator.hpp"

#include "hyteg/forms/P1RowSumForm.hpp"
#include "hyteg/forms/P1WrapperForm.hpp"
#include "hyteg/forms/form_fenics_base/P1ToP2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/forms/form_fenics_base/P2ToP1FenicsForm.hpp"
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

namespace hyteg {

using walberla::real_t;

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
: P1Operator< P1Form, Diagonal, Lumped, InvertDiagonal >( storage, minLevel, maxLevel, form )
{
   // pre-assemble edge, face and cell stencils
   assembleStencils();
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::scale( real_t scalar )
{
   // todo: Remove these constraints
   WALBERLA_CHECK_GREATER_EQUAL( minLevel_, 2, "scale() not implemented for level < 2" )
   WALBERLA_CHECK( !storage_->hasGlobalCells(), "scale() not implemented for macro-cells" )

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
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                                                                               const P1Function< idx_t >&                  src,
                                                                               const P1Function< idx_t >&                  dst,
                                                                               size_t                                      level,
                                                                               DoFType flag ) const
{
   for ( auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         vertexdof::macrovertex::saveOperator(
             vertex, this->getVertexStencilID(), src.getVertexDataID(), dst.getVertexDataID(), mat, level );
      }
   }
   if ( level >= 1 )
   {
      for ( auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            vertexdof::macroedge::saveOperator(
                level, edge, *storage_, this->getEdgeStencilID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat );
         }
      }
   }
   if ( level >= 2 )
   {
      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            if ( storage_->hasGlobalCells() )
            {
               vertexdof::macroface::saveOperator3D(
                   level, face, *storage_, this->getFaceStencil3DID(), src.getFaceDataID(), dst.getFaceDataID(), mat );
            }
            else
            {
               vertexdof::macroface::saveOperator(
                   level, face, this->getFaceStencilID(), src.getFaceDataID(), dst.getFaceDataID(), mat );
            }
         }
      }

      for ( auto& it : storage_->getCells() )
      {
         Cell& cell = *it.second;

         const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
         if ( testFlag( cellBC, flag ) )
         {
            vertexdof::macrocell::saveOperator(
                level, cell, this->getCellStencilID(), src.getCellDataID(), dst.getCellDataID(), mat );
         }
      }
   }
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
inline void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::apply_face3D_generated(
    Face&                                                    face,
    const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcId,
    const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
    const uint_t&                                            level,
    UpdateType                                               update ) const
{
   if ( face.getNumNeighborCells() == 2 )
   {
      WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start( "Two-sided" ); }
   }
   else
   {
      WALBERLA_NON_OPENMP_SECTION() { this->timingTree_->start( "One-sided" ); }
   }

   auto&        opr_data    = face.getData( faceStencil3DID_ )->getData( level );
   auto         src_data    = face.getData( srcId )->getPointer( level );
   auto         dst_data    = face.getData( dstId )->getPointer( level );
   const uint_t offset_gl_0 = levelinfo::num_microvertices_per_face( level );

   auto neighborCell0 = storage_->getCell( face.neighborCells()[0] );

   auto neighbor_cell_0_local_vertex_id_0 = static_cast< int32_t >(
       neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at( 0 ) );
   auto neighbor_cell_0_local_vertex_id_1 = static_cast< int32_t >(
       neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at( 1 ) );
   auto neighbor_cell_0_local_vertex_id_2 = static_cast< int32_t >(
       neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at( 2 ) );

   if ( update == Replace )
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
      vertexdof::macroface::generated::apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add( dst_data,
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
      const uint_t offset_gl_1 =
          offset_gl_0 + levelinfo::num_microvertices_per_face_from_width( levelinfo::num_microvertices_per_edge( level ) - 1 );

      auto neighborCell1 = storage_->getCell( face.neighborCells()[1] );

      auto neighbor_cell_1_local_vertex_id_0 = static_cast< int32_t >(
          neighborCell1->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell1->getLocalFaceID( face.getID() ) ).at( 0 ) );
      auto neighbor_cell_1_local_vertex_id_1 = static_cast< int32_t >(
          neighborCell1->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell1->getLocalFaceID( face.getID() ) ).at( 1 ) );
      auto neighbor_cell_1_local_vertex_id_2 = static_cast< int32_t >(
          neighborCell1->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell1->getLocalFaceID( face.getID() ) ).at( 2 ) );

      vertexdof::macroface::generated::apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add( dst_data,
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

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
inline void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::apply_face_generated(
    Face&                                                    face,
    const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcId,
    const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
    const uint_t&                                            level,
    UpdateType                                               update ) const
{
   real_t* opr_data = face.getData( faceStencilID_ )->getPointer( level );
   real_t* src_data = face.getData( srcId )->getPointer( level );
   real_t* dst_data = face.getData( dstId )->getPointer( level );

   if ( update == hyteg::Replace )
   {
      vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_replace(
          dst_data, src_data, opr_data, static_cast< int32_t >( level ) );
   }
   else if ( update == hyteg::Add )
   {
      vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_add(
          dst_data, src_data, opr_data, static_cast< int32_t >( level ) );
   }
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
inline void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::apply_cell_generated(
    Cell&                                                    cell,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >& srcId,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >& dstId,
    const uint_t&                                            level,
    UpdateType                                               update ) const
{
   auto&   opr_data = cell.getData( cellStencilID_ )->getData( level );
   real_t* src_data = cell.getData( srcId )->getPointer( level );
   real_t* dst_data = cell.getData( dstId )->getPointer( level );

   if ( update == Replace )
   {
      vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_replace(
          dst_data, src_data, static_cast< int32_t >( level ), opr_data );
   }
   else if ( update == Add )
   {
      vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_add(
          dst_data, src_data, static_cast< int32_t >( level ), opr_data );
   }
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
inline void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::smooth_sor_face3D_generated(
    Face&                                                    face,
    const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsId,
    const uint_t&                                            level,
    real_t                                                   relax,
    const bool&                                              backwards ) const
{
   auto  rhs_data = face.getData( rhsId )->getPointer( level );
   auto  dst_data = face.getData( dstId )->getPointer( level );
   auto& stencil  = face.getData( faceStencil3DID_ )->getData( level );

   auto neighborCell0 = storage_->getCell( face.neighborCells()[0] );

   auto neighbor_cell_0_local_vertex_id_0 = static_cast< int32_t >(
       neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at( 0 ) );
   auto neighbor_cell_0_local_vertex_id_1 = static_cast< int32_t >(
       neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at( 1 ) );
   auto neighbor_cell_0_local_vertex_id_2 = static_cast< int32_t >(
       neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at( 2 ) );

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

      auto neighbor_cell_1_local_vertex_id_0 = static_cast< int32_t >(
          neighborCell1->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell1->getLocalFaceID( face.getID() ) ).at( 0 ) );
      auto neighbor_cell_1_local_vertex_id_1 = static_cast< int32_t >(
          neighborCell1->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell1->getLocalFaceID( face.getID() ) ).at( 1 ) );
      auto neighbor_cell_1_local_vertex_id_2 = static_cast< int32_t >(
          neighborCell1->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell1->getLocalFaceID( face.getID() ) ).at( 2 ) );
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

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
inline void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::smooth_sor_face_generated(
    Face&                                                    face,
    const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsId,
    const uint_t&                                            level,
    real_t                                                   relax,
    const bool&                                              backwards ) const
{
   auto rhs_data = face.getData( rhsId )->getPointer( level );
   auto dst_data = face.getData( dstId )->getPointer( level );
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

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
inline void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::smooth_sor_cell_generated(
    Cell&                                                    cell,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >& dstId,
    const PrimitiveDataID< FunctionMemory< real_t >, Cell >& rhsId,
    const uint_t&                                            level,
    real_t                                                   relax,
    const bool&                                              backwards ) const
{
   auto  rhs_data = cell.getData( rhsId )->getPointer( level );
   auto  dst_data = cell.getData( dstId )->getPointer( level );
   auto& stencil  = cell.getData( cellStencilID_ )->getData( level );

   if ( backwards )
   {
      vertexdof::macrocell::generated::sor_3D_macrocell_P1_backwards(
          dst_data, rhs_data, static_cast< int32_t >( level ), stencil, relax );
   }
   else
   {
      vertexdof::macrocell::generated::sor_3D_macrocell_P1( dst_data, rhs_data, static_cast< int32_t >( level ), stencil, relax );
   }
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::assembleStencils()
{
   for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
   {
      if ( level >= 2 )
      {
         for ( const auto& it : storage_->getCells() )
         {
            auto  cell          = it.second;
            auto& stencilMemory = cell->getData( cellStencilID_ )->getData( level );
            assemble_variableStencil_cell_init( *cell, level );
            assemble_variableStencil_cell( stencilMemory, 1, 1, 1 );
         }

         for ( auto& it : storage_->getFaces() )
         {
            Face& face          = *it.second;
            auto  stencilMemory = face.getData( faceStencilID_ )->getPointer( level );
            auto& stencilMap    = face.getData( faceStencil3DID_ )->getData( level );

            assemble_variableStencil_face_init( face, level );

            if ( storage_->hasGlobalCells() )
            {
               assemble_variableStencil_face3D( stencilMap, 1, 1 );
               WALBERLA_ASSERT_GREATER( face.getNumNeighborCells(), 0 );
            }
            else
            {
               assemble_variableStencil_face( stencilMemory, 0, 0 );
            }
         }
      }
      if ( level >= 1 )
      {
         for ( auto& it : storage_->getEdges() )
         {
            Edge& edge          = *it.second;
            auto  stencilMemory = edge.getData( edgeStencilID_ )->getPointer( level );
            auto& stencilMap    = edge.getData( edgeStencil3DID_ )->getData( level );

            assemble_variableStencil_edge_init( edge, level );

            if ( storage_->hasGlobalCells() )
            {
               assemble_variableStencil_edge3D( stencilMap, 1 );
            }

            assemble_variableStencil_edge( stencilMemory, 1 ); // assemble both new and old version of 3D stencil
         }
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

template class P1ConstantOperator<
    P1WrapperForm< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > > >;
template class P1ConstantOperator<
    P1WrapperForm< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > > >;

template class P1ConstantOperator<
    P1WrapperForm< P2FenicsForm< p2_divt_cell_integral_0_otherwise, p2_tet_divt_tet_cell_integral_0_otherwise > > >;
template class P1ConstantOperator<
    P1WrapperForm< P2FenicsForm< p2_divt_cell_integral_1_otherwise, p2_tet_divt_tet_cell_integral_1_otherwise > > >;
template class P1ConstantOperator<
    P1WrapperForm< P2FenicsForm< fenics::NoAssemble, p2_tet_divt_tet_cell_integral_2_otherwise > > >;
template class P1ConstantOperator<
    P1WrapperForm< P2FenicsForm< p2_div_cell_integral_0_otherwise, p2_tet_div_tet_cell_integral_0_otherwise > > >;
template class P1ConstantOperator<
    P1WrapperForm< P2FenicsForm< p2_div_cell_integral_1_otherwise, p2_tet_div_tet_cell_integral_1_otherwise > > >;
template class P1ConstantOperator<
    P1WrapperForm< P2FenicsForm< fenics::NoAssemble, p2_tet_div_tet_cell_integral_2_otherwise > > >;

template class P1ConstantOperator< P1WrapperForm<
    P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_0_otherwise, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > > >;
template class P1ConstantOperator< P1WrapperForm<
    P1ToP2FenicsForm< p1_to_p2_divt_cell_integral_1_otherwise, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > > >;
template class P1ConstantOperator<
    P1WrapperForm< P1ToP2FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > > >;

template class P1ConstantOperator<
    P1WrapperForm< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > > >;
template class P1ConstantOperator<
    P1WrapperForm< P2ToP1FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > > >;
template class P1ConstantOperator<
    P1WrapperForm< P2ToP1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > > >;

template class P1ConstantOperator<
    P1WrapperForm< P2FenicsForm< p2_pspg_cell_integral_0_otherwise, p2_tet_pspg_tet_cell_integral_0_otherwise > > >;

template class P1ConstantOperator< P1LinearCombinationForm >;
template class P1ConstantOperator< P1WrapperForm< P2LinearCombinationForm > >;
template class P1ConstantOperator< P1WrapperForm< P2RowSumForm > >;

// Mostly for testing, as the P1ConstantOperator inherently was designed to support lumping
template class P1ConstantOperator< P1RowSumForm >;


// The following instantiations are required as building blocks in the P1ConstantEpsilonOperator class
// clang-format off
template class P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_0_otherwise, p1_tet_stokes_epsilon_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_1_otherwise, p1_tet_stokes_epsilon_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble                         , p1_tet_stokes_epsilon_tet_cell_integral_2_otherwise > >;

template class P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_2_otherwise, p1_tet_stokes_epsilon_tet_cell_integral_3_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_3_otherwise, p1_tet_stokes_epsilon_tet_cell_integral_4_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble                         , p1_tet_stokes_epsilon_tet_cell_integral_5_otherwise > >;

template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble                         , p1_tet_stokes_epsilon_tet_cell_integral_6_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble                         , p1_tet_stokes_epsilon_tet_cell_integral_7_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble                         , p1_tet_stokes_epsilon_tet_cell_integral_8_otherwise > >;
// clang-format on

// The following instantiations are required as building blocks in the P2ConstantEpsilon Operator_old class
// clang-format off
template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< p2_stokes_epsilon_cell_integral_0_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_0_otherwise > > >;
template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< p2_stokes_epsilon_cell_integral_1_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_1_otherwise > > >;
template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_2_otherwise > > >;

template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< p2_stokes_epsilon_cell_integral_2_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_3_otherwise > > >;
template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< p2_stokes_epsilon_cell_integral_3_otherwise, p2_tet_stokes_epsilon_tet_cell_integral_4_otherwise > > >;
template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_5_otherwise > > >;

template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_6_otherwise > > >;
template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_7_otherwise > > >;
template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< fenics::NoAssemble                         , p2_tet_stokes_epsilon_tet_cell_integral_8_otherwise > > >;
// clang-format on

// The following instantiations are required as building blocks in the P2ConstantFullViscousOperator_old class
// clang-format off
template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< p2_stokes_full_cell_integral_0_otherwise, p2_tet_stokes_full_tet_cell_integral_0_otherwise > > >;
template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< p2_stokes_full_cell_integral_1_otherwise, p2_tet_stokes_full_tet_cell_integral_1_otherwise > > >;
template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_2_otherwise > > >;

template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< p2_stokes_full_cell_integral_2_otherwise, p2_tet_stokes_full_tet_cell_integral_3_otherwise > > >;
template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< p2_stokes_full_cell_integral_3_otherwise, p2_tet_stokes_full_tet_cell_integral_4_otherwise > > >;
template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_5_otherwise > > >;

template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_6_otherwise > > >;
template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_7_otherwise > > >;
template class P1ConstantOperator< P1WrapperForm<P2FenicsForm< fenics::NoAssemble                      , p2_tet_stokes_full_tet_cell_integral_8_otherwise > > >;
// clang-format on
} // namespace hyteg
