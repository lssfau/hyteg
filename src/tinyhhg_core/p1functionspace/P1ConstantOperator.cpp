#include "P1ConstantOperator.hpp"

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/p1functionspace/generated/p1_diffusion.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include "tinyhhg_core/p1functionspace/VertexDoFMacroCell.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp"
#include "tinyhhg_core/p1functionspace/variablestencil/VertexDoFVariableStencil.hpp"

#include "tinyhhg_core/p2functionspace/generated_new/P2FenicsForm.hpp"

#include "P1Elements.hpp"
#include "generatedKernels/GeneratedKernelsVertexToVertexMacroCell3D.hpp"
#include "generatedKernels/GeneratedKernelsVertexToVertexMacroFace2D.hpp"

namespace hhg {

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::P1ConstantOperator(
    const std::shared_ptr< PrimitiveStorage >& storage,
    size_t                                     minLevel,
    size_t                                     maxLevel )
: Operator( storage, minLevel, maxLevel )
{
   auto cellP1StencilMemoryDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Cell > >(
       minLevel_, maxLevel_, vertexDoFMacroCellStencilMemorySize );
   auto faceP1StencilMemoryDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Face > >(
       minLevel_, maxLevel_, vertexDoFMacroFaceStencilMemorySize );
   auto edgeP1StencilMemoryDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Edge > >(
       minLevel_, maxLevel_, vertexDoFMacroEdgeStencilMemorySize );
   auto vertexP1StencilMemoryDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Vertex > >(
       minLevel_, maxLevel_, vertexDoFMacroVertexStencilMemorySize );

   storage->addCellData( cellStencilID_, cellP1StencilMemoryDataHandling, "P1OperatorCellStencil" );
   storage->addFaceData( faceStencilID_, faceP1StencilMemoryDataHandling, "P1OperatorFaceStencil" );
   storage->addEdgeData( edgeStencilID_, edgeP1StencilMemoryDataHandling, "P1OperatorEdgeStencil" );
   storage->addVertexData( vertexStencilID_, vertexP1StencilMemoryDataHandling, "P1OperatorVertexStencil" );

   if ( storage_->hasGlobalCells() )
   {
      const bool assemblyDefined = form.assembly3DDefined();
      WALBERLA_CHECK( assemblyDefined, "Assembly undefined for 3D elements." );
      if ( form.assemble3D() )
      {
         assembleStencils3D();
      }
   }
   else
   {
      if ( form.assemble2D() )
      {
         const bool assemblyDefined = form.assembly2DDefined();
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

         form.geometryMap = vertex->getGeometryMap();
         auto stencil =
             P1Elements::P1Elements3D::assembleP1LocalStencil( storage_, *vertex, indexing::Index( 0, 0, 0 ), level, form );

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

      for ( const auto& it : storage_->getEdges() )
      {
         auto edge          = it.second;
         auto stencilSize   = edge->getData( getEdgeStencilID() )->getSize( level );
         auto stencilMemory = edge->getData( getEdgeStencilID() )->getPointer( level );

         form.geometryMap = edge->getGeometryMap();
         auto stencil =
             P1Elements::P1Elements3D::assembleP1LocalStencil( storage_, *edge, indexing::Index( 1, 0, 0 ), level, form );

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
               stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell, it.second->getNumNeighborFaces() )] =
                   0;
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
               stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell, it.second->getNumNeighborFaces() )] =
                   0;
            }
         }
         if ( InvertDiagonal )
         {
            stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )] =
                1.0 / stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )];
         }
      }

      for ( const auto& it : storage_->getFaces() )
      {
         auto face          = it.second;
         auto stencilMemory = face->getData( getFaceStencilID() )->getPointer( level );

         form.geometryMap = face->getGeometryMap();
         auto stencil =
             P1Elements::P1Elements3D::assembleP1LocalStencil( storage_, *face, indexing::Index( 1, 1, 0 ), level, form );

         if ( face->getNumNeighborCells() == 1 )
         {
            for ( const auto stencilDir : vertexdof::macroface::neighborsWithOneNeighborCellWithCenter )
            {
               if ( stencil.count( stencilDir ) == 0 )
               {
                  stencil[stencilDir] = real_c( 0 );
               }
            }
         }
         else
         {
            for ( const auto stencilDir : vertexdof::macroface::neighborsWithTwoNeighborCellsWithCenter )
            {
               if ( stencil.count( stencilDir ) == 0 )
               {
                  stencil[stencilDir] = real_c( 0 );
               }
            }
         }

         for ( const auto stencilIt : stencil )
         {
            const auto stencilIdx     = vertexdof::stencilIndexFromVertex( stencilIt.first );
            stencilMemory[stencilIdx] = stencil[stencilIt.first];
         }

         if ( Lumped )
         {
            if ( face->getNumNeighborCells() == 1 )
            {
               for ( auto dir : vertexdof::macroface::neighborsWithOneNeighborCellWithoutCenter )
               {
                  stencilMemory[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] +=
                      stencilMemory[vertexdof::stencilIndexFromVertex( dir )];
                  stencilMemory[vertexdof::stencilIndexFromVertex( dir )] = 0;
               }
            }
            else
            {
               for ( auto dir : vertexdof::macroface::neighborsWithTwoNeighborCellsWithoutCenter )
               {
                  stencilMemory[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] +=
                      stencilMemory[vertexdof::stencilIndexFromVertex( dir )];
                  stencilMemory[vertexdof::stencilIndexFromVertex( dir )] = 0;
               }
            }
         }
         if ( Diagonal )
         {
            if ( face->getNumNeighborCells() == 1 )
            {
               for ( auto dir : vertexdof::macroface::neighborsWithOneNeighborCellWithoutCenter )
               {
                  stencilMemory[vertexdof::stencilIndexFromVertex( dir )] = 0;
               }
            }
            else
            {
               for ( auto dir : vertexdof::macroface::neighborsWithTwoNeighborCellsWithoutCenter )
               {
                  stencilMemory[vertexdof::stencilIndexFromVertex( dir )] = 0;
               }
            }
         }
         if ( InvertDiagonal )
         {
            stencilMemory[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] =
                1.0 / stencilMemory[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];
         }
      }

      for ( const auto& it : storage_->getCells() )
      {
         auto cell          = it.second;
         auto stencilMemory = cell->getData( getCellStencilID() )->getPointer( level );

         form.geometryMap = cell->getGeometryMap();
         auto stencil =
             P1Elements::P1Elements3D::assembleP1LocalStencil( storage_, *cell, indexing::Index( 1, 1, 1 ), level, form );

         for ( const auto stencilIt : stencil )
         {
            const auto stencilIdx     = vertexdof::stencilIndexFromVertex( stencilIt.first );
            stencilMemory[stencilIdx] = stencil[stencilIt.first];
         }

         if ( Lumped )
         {
            for ( auto dir : vertexdof::macrocell::neighborsWithoutCenter )
            {
               stencilMemory[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] +=
                   stencilMemory[vertexdof::stencilIndexFromVertex( dir )];
               stencilMemory[vertexdof::stencilIndexFromVertex( dir )] = 0;
            }
         }
         if ( Diagonal )
         {
            for ( auto dir : vertexdof::macrocell::neighborsWithoutCenter )
            {
               stencilMemory[vertexdof::stencilIndexFromVertex( dir )] = 0;
            }
         }
         if ( InvertDiagonal )
         {
            stencilMemory[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] =
                1.0 / stencilMemory[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];
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

         form.geometryMap = face.getGeometryMap();

         Point3D dirS  = -1.0 * d2;
         Point3D dirSE = d0 - 1.0 * d2;
         Point3D dirE  = d0;
         Point3D dirW  = -1.0 * d0;
         Point3D dirNW = -1.0 * d0 + d2;
         Point3D dirN  = d2;

         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form, {x, x + dirW, x + dirS}, P1Elements::P1Elements2D::elementSW, face_stencil );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form, {x, x + dirS, x + dirSE}, P1Elements::P1Elements2D::elementS, face_stencil );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form, {x, x + dirSE, x + dirE}, P1Elements::P1Elements2D::elementSE, face_stencil );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form, {x, x + dirE, x + dirN}, P1Elements::P1Elements2D::elementNE, face_stencil );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form, {x, x + dirN, x + dirNW}, P1Elements::P1Elements2D::elementN, face_stencil );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form, {x, x + dirNW, x + dirW}, P1Elements::P1Elements2D::elementNW, face_stencil );

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
         form.geometryMap = faceS->getGeometryMap();
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form, {x, x + dir_W, x + dir_S}, P1Elements::P1Elements2D::elementSW, edge_stencil );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form, {x, x + dir_S, x + dir_SE}, P1Elements::P1Elements2D::elementS, edge_stencil );
         vertexdof::variablestencil::assembleLocalStencil< P1Form >(
             form, {x, x + dir_SE, x + dir_E}, P1Elements::P1Elements2D::elementSE, edge_stencil );

         if ( edge.getNumNeighborFaces() == 2 )
         {
            form.geometryMap = faceN->getGeometryMap();
            vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                form, {x, x + dir_E, x + dir_N}, P1Elements::P1Elements2D::elementNE, edge_stencil );
            vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                form, {x, x + dir_N, x + dir_NW}, P1Elements::P1Elements2D::elementN, edge_stencil );
            vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                form, {x, x + dir_NW, x + dir_W}, P1Elements::P1Elements2D::elementNW, edge_stencil );
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
            Face* face       = storage_->getFace( faceId );
            form.geometryMap = face->getGeometryMap();

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
            form.integrate( {{x, x + d0, x + d2}}, matrixRow );

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
   this->startTiming( "Apply" );
   src.communicate< Vertex, Edge >( level );
   src.communicate< Edge, Face >( level );
   src.communicate< Face, Cell >( level );

   src.communicate< Cell, Face >( level );
   src.communicate< Face, Edge >( level );
   src.communicate< Edge, Vertex >( level );

   for ( const auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         vertexdof::macrovertex::apply< real_t >(
             vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), level, updateType );
      }
   }

   for ( const auto& it : storage_->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         vertexdof::macroedge::apply< real_t >(
             level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
      }
   }

   for ( const auto& it : storage_->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         if ( hhg::globalDefines::useGeneratedKernels && ( !storage_->hasGlobalCells() ) )
         {
            real_t* opr_data = face.getData( faceStencilID_ )->getPointer( level );
            real_t* src_data = face.getData( src.getFaceDataID() )->getPointer( level );
            real_t* dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );
            if ( updateType == hhg::Replace )
            {
               vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_replace(
                   dst_data, src_data, opr_data, static_cast< int64_t >( level ) );
            }
            else if ( updateType == hhg::Add )
            {
               vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_add(
                   dst_data, src_data, opr_data, static_cast< int64_t >( level ) );
            }
         }
         else
         {
            vertexdof::macroface::apply< real_t >(
                level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
         }
      }
   }

   for ( const auto& it : storage_->getCells() )
   {
      Cell& cell = *it.second;

      const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
      if ( testFlag( cellBC, flag ) )
      {
         if ( hhg::globalDefines::useGeneratedKernels )
         {
            real_t* opr_data = cell.getData( cellStencilID_ )->getPointer( level );
            real_t* src_data = cell.getData( src.getCellDataID() )->getPointer( level );
            real_t* dst_data = cell.getData( dst.getCellDataID() )->getPointer( level );
            if ( updateType == Replace )
            {
               vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_replace(
                   dst_data, src_data, opr_data, static_cast< int64_t >( level ) );
            }
            else if ( updateType == Add )
            {
               vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_add(
                   dst_data, src_data, opr_data, static_cast< int64_t >( level ) );
            }
         }
         else
         {
            vertexdof::macrocell::apply< real_t >(
                level, cell, cellStencilID_, src.getCellDataID(), dst.getCellDataID(), updateType );
         }
      }
   }
   this->stopTiming( "Apply" );
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::smooth_gs( const P1Function< real_t >& dst,
                                                                                const P1Function< real_t >& rhs,
                                                                                size_t                      level,
                                                                                DoFType                     flag ) const
{
   dst.communicate< Vertex, Edge >( level );
   dst.communicate< Edge, Face >( level );
   dst.communicate< Face, Cell >( level );

   dst.communicate< Cell, Face >( level );
   dst.communicate< Face, Edge >( level );
   dst.communicate< Edge, Vertex >( level );

   for ( auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         vertexdof::macrovertex::smooth_gs< real_t >(
             vertex, vertexStencilID_, dst.getVertexDataID(), rhs.getVertexDataID(), level );
      }
   }

   dst.communicate< Vertex, Edge >( level );

   for ( auto& it : storage_->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         vertexdof::macroedge::smooth_gs< real_t >( level, edge, edgeStencilID_, dst.getEdgeDataID(), rhs.getEdgeDataID() );
      }
   }

   dst.communicate< Edge, Face >( level );

   for ( auto& it : storage_->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         if ( hhg::globalDefines::useGeneratedKernels && ( !storage_->hasGlobalCells() ) )
         {
            real_t* opr_data = face.getData( faceStencilID_ )->getPointer( level );
            real_t* dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );
            real_t* rhs_data = face.getData( rhs.getFaceDataID() )->getPointer( level );
            vertexdof::macroface::generated::gaussseidel_2D_macroface_vertexdof_to_vertexdof(
                dst_data, rhs_data, opr_data, static_cast< int64_t >( level ) );
         }
         else
         {
            vertexdof::macroface::smooth_gs< real_t >( level, face, faceStencilID_, dst.getFaceDataID(), rhs.getFaceDataID() );
         }
      }
   }

   dst.communicate< Face, Cell >( level );

   for ( auto& it : storage_->getCells() )
   {
      Cell& cell = *it.second;

      const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
      if ( testFlag( cellBC, flag ) )
      {
         vertexdof::macrocell::smooth_gs< real_t >( level, cell, cellStencilID_, dst.getCellDataID(), rhs.getCellDataID() );
      }
   }
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::smooth_sor( const P1Function< real_t >& dst,
                                                                                 const P1Function< real_t >& rhs,
                                                                                 real_t                      relax,
                                                                                 size_t                      level,
                                                                                 DoFType                     flag ) const
{
   dst.communicate< Vertex, Edge >( level );
   dst.communicate< Edge, Face >( level );
   dst.communicate< Face, Cell >( level );

   dst.communicate< Cell, Face >( level );
   dst.communicate< Face, Edge >( level );
   dst.communicate< Edge, Vertex >( level );

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

   dst.communicate< Vertex, Edge >( level );

   for ( auto& it : storage_->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         vertexdof::macroedge::smooth_sor< real_t >(
             level, edge, edgeStencilID_, dst.getEdgeDataID(), rhs.getEdgeDataID(), relax );
      }
   }

   dst.communicate< Edge, Face >( level );

   for ( auto& it : storage_->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         vertexdof::macroface::smooth_sor< real_t >(
             level, face, faceStencilID_, dst.getFaceDataID(), rhs.getFaceDataID(), relax );
      }
   }

   dst.communicate< Face, Cell >( level );

   for ( auto& it : storage_->getCells() )
   {
      Cell& cell = *it.second;

      const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
      if ( testFlag( cellBC, flag ) )
      {
         vertexdof::macrocell::smooth_sor< real_t >(
             level, cell, cellStencilID_, dst.getCellDataID(), rhs.getCellDataID(), relax );
      }
   }
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::smooth_jac( const P1Function< real_t >& dst,
                                                                                 const P1Function< real_t >& rhs,
                                                                                 const P1Function< real_t >& tmp,
                                                                                 size_t                      level,
                                                                                 DoFType                     flag ) const
{
   tmp.communicate< Vertex, Edge >( level );
   tmp.communicate< Edge, Face >( level );
   tmp.communicate< Face, Cell >( level );

   tmp.communicate< Cell, Face >( level );
   tmp.communicate< Face, Edge >( level );
   tmp.communicate< Edge, Vertex >( level );

   for ( auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if ( testFlag( vertexBC, flag ) )
      {
         vertexdof::macrovertex::smooth_jac(
             vertex, vertexStencilID_, dst.getVertexDataID(), rhs.getVertexDataID(), tmp.getVertexDataID(), level );
      }
   }

   for ( auto& it : storage_->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if ( testFlag( edgeBC, flag ) )
      {
         vertexdof::macroedge::smooth_jac< real_t >(
             level, edge, edgeStencilID_, dst.getEdgeDataID(), rhs.getEdgeDataID(), tmp.getEdgeDataID() );
      }
   }

   for ( auto& it : storage_->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if ( testFlag( faceBC, flag ) )
      {
         vertexdof::macroface::smooth_jac< real_t >(
             level, face, faceStencilID_, dst.getFaceDataID(), rhs.getFaceDataID(), tmp.getFaceDataID() );
      }
   }
}

template < class P1Form, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< P1Form, Diagonal, Lumped, InvertDiagonal >::scale( real_t scalar )
{
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

template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, fenics::NoAssemble > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, fenics::UndefinedAssembly> >;

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
                                   true >;

template class P1ConstantOperator< P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise > >;

template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::UndefinedAssembly, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >;

template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > >;
template class P1ConstantOperator< P1FenicsForm< fenics::UndefinedAssembly, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > >;

template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_tet_diffusion_cell_integral_0_otherwise> >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_tet_mass_cell_integral_0_otherwise> >;

template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_divt_cell_integral_0_otherwise> >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_divt_cell_integral_1_otherwise> >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_div_cell_integral_0_otherwise> >;
template class P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_div_cell_integral_1_otherwise> >;

template class P1ConstantOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >;

template class P1ConstantOperator< P2FenicsForm< p2_divt_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< p2_divt_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< p2_div_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< p2_div_cell_integral_1_otherwise > >;

template class P1ConstantOperator< P2FenicsForm< p1_to_p2_divt_cell_integral_0_otherwise, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< p1_to_p2_divt_cell_integral_1_otherwise, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< fenics::NoAssemble,                      p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > >;

template class P1ConstantOperator< P2FenicsForm< p2_to_p1_div_cell_integral_0_otherwise, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< p2_to_p1_div_cell_integral_1_otherwise, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > >;
template class P1ConstantOperator< P2FenicsForm< fenics::NoAssemble,                     p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >;

} // namespace hhg
