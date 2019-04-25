#include "P1ConstantOperator.hpp"

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/p1functionspace/generated/p1_diffusion.h"
#include "tinyhhg_core/p1functionspace/generated/p1_div.h"
#include "tinyhhg_core/p1functionspace/generated/p1_divt.h"
#include "tinyhhg_core/p1functionspace/generated/p1_mass.h"
#include "tinyhhg_core/p1functionspace/generated/p1_pspg.h"
#include "tinyhhg_core/p1functionspace/generated/p1_stokes_epsilon.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_diffusion.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_div_tet.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_divt_tet.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_mass.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_pspg_tet.h"

#include "tinyhhg_core/p2functionspace/generated/p2_tet_diffusion.h"
#include "tinyhhg_core/p2functionspace/generated/p2_tet_mass.h"
#include "tinyhhg_core/p2functionspace/generated/p2_tet_div_tet.h"
#include "tinyhhg_core/p2functionspace/generated/p2_tet_divt_tet.h"
#include "tinyhhg_core/p2functionspace/generated/p2_tet_pspg_tet.h"

#include "tinyhhg_core/mixedoperators/generated/p2_to_p1_div.h"
#include "tinyhhg_core/mixedoperators/generated/p2_to_p1_tet_div_tet.h"
#include "tinyhhg_core/mixedoperators/generated/p1_to_p2_tet_divt_tet.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include "generatedKernels/GeneratedKernelsVertexToVertexMacroFace2D.hpp"
#include "generatedKernels/GeneratedKernelsVertexToVertexMacroCell3D.hpp"
#include "generatedKernels/GeneratedKernelsVertexToVertexMacroFace3D.hpp"
#include "P1Elements.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroCell.hpp"
#include "tinyhhg_core/LevelWiseMemory.hpp"

namespace hhg {

template < class UFCOperator2D, class UFCOperator3D, bool Diagonal, bool Lumped, bool InvertDiagonal >
P1ConstantOperator< UFCOperator2D, UFCOperator3D, Diagonal, Lumped, InvertDiagonal >::P1ConstantOperator(
    const std::shared_ptr< PrimitiveStorage >& storage,
    size_t                                     minLevel,
    size_t                                     maxLevel )
: Operator( storage, minLevel, maxLevel )
{
   auto cellP1StencilMemoryDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell > >(
           minLevel_, maxLevel_ );

   auto face3DP1StencilMemoryDataHandling =
       std::make_shared< LevelWiseMemoryDataHandling< LevelWiseMemory< vertexdof::macroface::StencilMap_T >, Face > >(
           minLevel_, maxLevel_ );

   auto faceP1StencilMemoryDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Face > >(
       minLevel_, maxLevel_, vertexDoFMacroFaceStencilMemorySize );
   auto edgeP1StencilMemoryDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Edge > >(
       minLevel_, maxLevel_, vertexDoFMacroEdgeStencilMemorySize );
   auto vertexP1StencilMemoryDataHandling = std::make_shared< MemoryDataHandling< StencilMemory< real_t >, Vertex > >(
       minLevel_, maxLevel_, vertexDoFMacroVertexStencilMemorySize );

   storage->addCellData( cellStencilID_, cellP1StencilMemoryDataHandling, "P1OperatorCellStencil" );
   storage->addFaceData( faceStencilID_, faceP1StencilMemoryDataHandling, "P1OperatorFaceStencil" );
   storage->addFaceData( faceStencil3DID_, face3DP1StencilMemoryDataHandling, "P1OperatorFaceStencil" );
   storage->addEdgeData( edgeStencilID_, edgeP1StencilMemoryDataHandling, "P1OperatorEdgeStencil" );
   storage->addVertexData( vertexStencilID_, vertexP1StencilMemoryDataHandling, "P1OperatorVertexStencil" );

   if ( storage_->hasGlobalCells() )
   {
      const bool assemblyDefined = !std::is_same< UFCOperator3D, hhg::fenics::UndefinedAssembly >::value;
      WALBERLA_CHECK( assemblyDefined, "Assembly undefined for 3D elements." );
      if ( !std::is_same< UFCOperator3D, fenics::NoAssemble >::value )
      {
         assembleStencils3D();
      }
   }
   else
   {
      if ( !std::is_same< UFCOperator2D, fenics::NoAssemble >::value )
      {
         const bool assemblyDefined = !std::is_same< UFCOperator2D, hhg::fenics::UndefinedAssembly >::value;
         WALBERLA_CHECK( assemblyDefined, "Assembly undefined for 2D elements." );
         assembleStencils();
      }
   }
}

template < class UFCOperator2D, class UFCOperator3D, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< UFCOperator2D, UFCOperator3D, Diagonal, Lumped, InvertDiagonal >::assembleStencils3D()
{
   for( uint_t level = minLevel_; level <= maxLevel_; level++ )
   {
      for( const auto& it : storage_->getVertices() )
      {
         auto          vertex        = it.second;
         auto          stencilSize   = vertex->getData( getVertexStencilID() )->getSize( level );
         auto          stencilMemory = vertex->getData( getVertexStencilID() )->getPointer( level );
         UFCOperator3D ufcOperator;

         auto stencil = P1Elements::P1Elements3D::assembleP1LocalStencil(
             storage_, *vertex, indexing::Index( 0, 0, 0 ), level, ufcOperator );

         WALBERLA_ASSERT_EQUAL( stencilSize, stencil.size() );
         for( uint_t i = 0; i < stencilSize; i++ )
         {
            stencilMemory[i] = stencil[i];
         }

         if( Lumped )
         {
            for( uint_t i = 1; i < stencilSize; i++ )
            {
               stencilMemory[0] += stencilMemory[i];
               stencilMemory[i] = 0;
            }
         }
         if( InvertDiagonal )
         {
            stencilMemory[0] = 1.0 / stencilMemory[0];
         }
      }

      for( const auto& it : storage_->getEdges() )
      {
         auto          edge          = it.second;
         auto          stencilSize   = edge->getData( getEdgeStencilID() )->getSize( level );
         auto          stencilMemory = edge->getData( getEdgeStencilID() )->getPointer( level );
         UFCOperator3D ufcOperator;

         auto stencil =
             P1Elements::P1Elements3D::assembleP1LocalStencil( storage_, *edge, indexing::Index( 1, 0, 0 ), level, ufcOperator );

         WALBERLA_ASSERT_EQUAL( stencilSize, stencil.size() );
         for( uint_t i = 0; i < stencilSize; i++ )
         {
            stencilMemory[i] = stencil[i];
         }
         if( Lumped )
         {
            stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )] +=
                stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_W )];
            stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_W )] = 0;
            stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )] +=
                stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_E )];
            stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_E )] = 0;
            for( uint_t neighborFace = 0; neighborFace < it.second->getNumNeighborFaces(); neighborFace++ )
            {
               stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )] +=
                   stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborFace( stencilDirection::VERTEX_W, neighborFace )];
               stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborFace( stencilDirection::VERTEX_W, neighborFace )] = 0;
               stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )] +=
                   stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborFace( stencilDirection::VERTEX_E, neighborFace )];
               stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborFace( stencilDirection::VERTEX_E, neighborFace )] = 0;
            }
            for( uint_t neighborCell = 0; neighborCell < it.second->getNumNeighborCells(); neighborCell++ )
            {
               stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )] +=
                   stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell,
                                                                                   it.second->getNumNeighborFaces() )];
               stencilMemory[vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell, it.second->getNumNeighborFaces() )] =
                   0;
            }
         }
         if( InvertDiagonal )
         {
            stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )] =
                1.0 / stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )];
         }
      }

      for ( const auto& it : storage_->getFaces() )
      {
         auto          face          = it.second;
         auto&         stencilMemory = face->getData( getFaceStencil3DID() )->getData( level );
         UFCOperator3D ufcOperator;

         for ( uint_t neighborCellID = 0; neighborCellID < face->getNumNeighborCells(); neighborCellID++ )
         {
            auto neighborCell = storage_->getCell( face->neighborCells().at( neighborCellID ) );
            auto vertexAssemblyIndexInCell =
                vertexdof::macroface::getIndexInNeighboringMacroCell( {1, 1, 0}, *face, neighborCellID, *storage_, level );
            stencilMemory[neighborCellID] = P1Elements::P1Elements3D::assembleP1LocalStencilNew(
                storage_, *neighborCell, vertexAssemblyIndexInCell, level, ufcOperator );
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
                  if ( !( neighborCellID == 0 && stencilIt.first == indexing::IndexIncrement( {0, 0, 0} ) ) )
                  {
                     stencilMemory[0][{0, 0, 0}] += stencilIt.second;
                     stencilIt.second = 0;
                  }
               }
            }
         }

         if ( InvertDiagonal )
         {
            for ( uint_t neighborCellID = 1; neighborCellID < face->getNumNeighborCells(); neighborCellID++ )
            {
               stencilMemory[0][{0, 0, 0}] += stencilMemory[neighborCellID][{0, 0, 0}];
               stencilMemory[neighborCellID][{0, 0, 0}] = 0;
            }
            stencilMemory[0][{0, 0, 0}] = 1.0 / stencilMemory[0][{0, 0, 0}];
         }
      }

      for( const auto& it : storage_->getCells() )
      {
         auto          cell          = it.second;
         auto &        stencilMemory = cell->getData( getCellStencilID() )->getData( level );
         UFCOperator3D ufcOperator;

         stencilMemory =
             P1Elements::P1Elements3D::assembleP1LocalStencilNew( storage_, *cell, indexing::Index( 1, 1, 1 ), level, ufcOperator );

         if( Lumped )
         {
            for( auto dir : vertexdof::macrocell::neighborsWithoutCenter )
            {
               stencilMemory[{ 0, 0, 0 }] +=
                   stencilMemory[vertexdof::logicalIndexOffsetFromVertex( dir )];
               stencilMemory[vertexdof::logicalIndexOffsetFromVertex( dir )] = 0;
            }
         }

         if( InvertDiagonal )
         {
            stencilMemory[{ 0, 0, 0 }] =
                1.0 / stencilMemory[{ 0, 0, 0 }];
         }
      }
   }
}

template < class UFCOperator2D, class UFCOperator3D, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< UFCOperator2D, UFCOperator3D, Diagonal, Lumped, InvertDiagonal >::compute_local_stiffness(
    const Face&         face,
    size_t              level,
    Matrix3r&           local_stiffness,
    fenics::ElementType element_type )
{
   real_t coords[6];
   fenics::compute_micro_coords( face, level, coords, element_type );
   UFCOperator2D gen;
   gen.tabulate_tensor( local_stiffness.data(), nullptr, coords, 0 );

   if( Diagonal )
   {
      for( size_t i = 0; i < 3; ++i )
      {
         for( size_t j = 0; j < 3; ++j )
         {
            if( i != j )
            {
               local_stiffness( i, j ) = real_t( 0 );
            }
         }
      }
   }
}

template < class UFCOperator2D, class UFCOperator3D, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< UFCOperator2D, UFCOperator3D, Diagonal, Lumped, InvertDiagonal >::assembleStencils()
{
   using namespace P1Elements::P1Elements2D;
   typedef stencilDirection sD;

   Matrix3r local_stiffness_gray;
   Matrix3r local_stiffness_blue;

   for( uint_t level = minLevel_; level <= maxLevel_; ++level )
   {
      for( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         auto face_stencil = face.getData( faceStencilID_ )->getPointer( level );
         compute_local_stiffness( face, level, local_stiffness_gray, fenics::GRAY );
         compute_local_stiffness( face, level, local_stiffness_blue, fenics::BLUE );

         for( uint_t i = 0; i < P1GrayElements.size(); ++i )
         {
            assembleP1LocalStencil( P1GrayStencilMaps[i], P1GrayDoFMaps[i], local_stiffness_gray, face_stencil );
         }

         for( uint_t i = 0; i < P1BlueElements.size(); ++i )
         {
            assembleP1LocalStencil( P1BlueStencilMaps[i], P1BlueDoFMaps[i], local_stiffness_blue, face_stencil );
         }

         if( Lumped )
         {
            for( const auto& neighbor : vertexdof::macroface::neighborsWithoutCenter )
            {
               face_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] +=
                   face_stencil[vertexdof::stencilIndexFromVertex( neighbor )];
               face_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
            }
         }

         if( InvertDiagonal )
         {
            face_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] =
                1.0 / face_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )];
         }
      }

      for( auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         auto edge_stencil = edge.getData( edgeStencilID_ )->getPointer( level );

         Face* firstFace = storage_->getFace( edge.neighborFaces()[0] );
         compute_local_stiffness( *firstFace, level, local_stiffness_gray, fenics::GRAY );
         compute_local_stiffness( *firstFace, level, local_stiffness_blue, fenics::BLUE );

         size_t start_id    = firstFace->vertex_index( edge.neighborVertices()[0] );
         size_t end_id      = firstFace->vertex_index( edge.neighborVertices()[1] );
         size_t opposite_id = firstFace->vertex_index( firstFace->get_vertex_opposite_to_edge( edge.getID() ) );

         assembleP1LocalStencil( convertStencilDirectionsToIndices( elementSW ),
                                 {{end_id, start_id, opposite_id}},
                                 local_stiffness_gray,
                                 edge_stencil );
         assembleP1LocalStencil( convertStencilDirectionsToIndices( elementS ),
                                 {{opposite_id, end_id, start_id}},
                                 local_stiffness_blue,
                                 edge_stencil );
         assembleP1LocalStencil( convertStencilDirectionsToIndices( elementSE ),
                                 {{start_id, opposite_id, end_id}},
                                 local_stiffness_gray,
                                 edge_stencil );

         if( edge.getNumNeighborFaces() == 2 )
         {
            Face* secondFace = storage_->getFace( edge.neighborFaces()[1] );
            compute_local_stiffness( *secondFace, level, local_stiffness_gray, fenics::GRAY );
            compute_local_stiffness( *secondFace, level, local_stiffness_blue, fenics::BLUE );

            size_t startIdSecondFace    = secondFace->vertex_index( edge.neighborVertices()[0] );
            size_t endIdSecondFace      = secondFace->vertex_index( edge.neighborVertices()[1] );
            size_t oppositeIdSecondFace = secondFace->vertex_index( secondFace->get_vertex_opposite_to_edge( edge.getID() ) );

            assembleP1LocalStencil( convertStencilDirectionsToIndices( elementNE ),
                                    {{startIdSecondFace, endIdSecondFace, oppositeIdSecondFace}},
                                    local_stiffness_gray,
                                    edge_stencil );
            assembleP1LocalStencil( convertStencilDirectionsToIndices( elementN ),
                                    {{oppositeIdSecondFace, startIdSecondFace, endIdSecondFace}},
                                    local_stiffness_blue,
                                    edge_stencil );
            assembleP1LocalStencil( convertStencilDirectionsToIndices( elementNW ),
                                    {{endIdSecondFace, oppositeIdSecondFace, startIdSecondFace}},
                                    local_stiffness_gray,
                                    edge_stencil );
         }

         if( Lumped )
         {
            for( const auto& neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
            {
               edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] +=
                   edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )];
               edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
            }

            for( const auto& neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
            {
               edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] +=
                   edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )];
               edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
            }

            if( edge.getNumNeighborFaces() == 2 )
            {
               for( const auto& neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
               {
                  edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] +=
                      edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )];
                  edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
               }
            }
         }

         if( InvertDiagonal )
         {
            edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] =
                1.0 / edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )];
         }

         //        WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("EDGE.id = {}:edge_stencil = {}", edge.getID().getID(), PointND<real_t, 7>(&edge_stencil[0])));
      }

      for( auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         auto vertex_stencil = vertex.getData( vertexStencilID_ )->getPointer( level );

         // iterate over adjacent faces
         for( auto& faceId : vertex.neighborFaces() )
         {
            Face* face = storage_->getFace( faceId );

            compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );

            uint_t v_i = face->vertex_index( vertex.getID() );

            std::vector< PrimitiveID > adj_edges = face->adjacent_edges( vertex.getID() );

            std::array< uint_t, 3 > stencilMap{};
            stencilMap[0] = 0;

            std::array< uint_t, 3 > dofMap{};
            dofMap[0] = v_i;

            // iterate over adjacent edges
            for( uint_t i = 0; i < adj_edges.size(); ++i )
            {
               uint_t      edge_idx = vertex.edge_index( adj_edges[i] ) + 1;
               Edge*       edge     = storage_->getEdge( adj_edges[i] );
               PrimitiveID vertex_j = edge->get_opposite_vertex( vertex.getID() );

               stencilMap[i + 1] = edge_idx;
               dofMap[i + 1]     = face->vertex_index( vertex_j );
            }

            assembleP1LocalStencil( stencilMap, dofMap, local_stiffness_gray, vertex_stencil );

            //          WALBERLA_LOG_DEVEL_ON_ROOT(fmt::format("VERTEX.id = {}:vertex_stencil = {}", vertex.getID().getID(), PointND<real_t, 3>(&vertex_stencil[0])));
         }

         if( Lumped )
         {
            for( uint_t i = 1; i < vertex.getData( vertexStencilID_ )->getSize( level ); ++i )
            {
               vertex_stencil[0] += vertex_stencil[i];
               vertex_stencil[i] = 0;
            }
         }

         if( InvertDiagonal )
         {
            vertex_stencil[0] = 1.0 / vertex_stencil[0];
         }
      }
   }
}

template<class UFCOperator2D, class UFCOperator3D, bool Diagonal, bool Lumped, bool InvertDiagonal>
void P1ConstantOperator<UFCOperator2D, UFCOperator3D, Diagonal, Lumped, InvertDiagonal>::apply(const P1Function<real_t> &src,
                                                                                               const P1Function<real_t> &dst,
                                                                                               size_t level,
                                                                                               DoFType flag,
                                                                                               UpdateType updateType) const{
   this->startTiming( "Apply" );
   src.communicate< Vertex, Edge >( level );
   src.communicate< Edge, Face >( level );
   src.communicate< Face, Cell >( level );

   src.communicate< Cell, Face >( level );
   src.communicate< Face, Edge >( level );
   src.communicate< Edge, Vertex >( level );

  this->timingTree_->start( "Macro-Vertex" );

   for( const auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if( testFlag( vertexBC, flag ) )
      {
         vertexdof::macrovertex::apply< real_t >(
                 vertex, vertexStencilID_, src.getVertexDataID(), dst.getVertexDataID(), level, updateType );
      }
   }

  this->timingTree_->stop( "Macro-Vertex" );

  this->timingTree_->start( "Macro-Edge" );

   for( const auto& it : storage_->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if( testFlag( edgeBC, flag ) )
      {
         vertexdof::macroedge::apply< real_t >(
                 level, edge, edgeStencilID_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
      }
   }

  this->timingTree_->stop( "Macro-Edge" );

  this->timingTree_->start( "Macro-Face" );

  for ( const auto& it : storage_->getFaces() )
  {
     Face& face = *it.second;

     const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
     if ( testFlag( faceBC, flag ) )
     {
        if ( storage_->hasGlobalCells() )
        {
           if ( hhg::globalDefines::useGeneratedKernels && face.getNumNeighborCells() == 2 && updateType == Replace )
           {
             auto opr_data = face.getData( faceStencil3DID_ )->getData( level );
             auto src_data = face.getData( src.getFaceDataID() )->getPointer( level );
             auto dst_data = face.getData( dst.getFaceDataID() )->getPointer( level );
             const uint_t offset_gl_0 = levelinfo::num_microvertices_per_face( level );
             const uint_t offset_gl_1 = offset_gl_0 + levelinfo::num_microvertices_per_face_from_width( levelinfo::num_microvertices_per_edge(level) - 1 );

             auto neighborCell0 = storage_->getCell( face.neighborCells()[0] );
             auto neighborCell1 = storage_->getCell( face.neighborCells()[1] );

             auto neighbor_cell_0_local_vertex_id_0 = static_cast< int64_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at(0) );
             auto neighbor_cell_0_local_vertex_id_1 = static_cast< int64_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at(1) );
             auto neighbor_cell_0_local_vertex_id_2 = static_cast< int64_t >( neighborCell0->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell0->getLocalFaceID( face.getID() ) ).at(2) );

             auto neighbor_cell_1_local_vertex_id_0 = static_cast< int64_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell1->getLocalFaceID( face.getID() ) ).at(0) );
             auto neighbor_cell_1_local_vertex_id_1 = static_cast< int64_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell1->getLocalFaceID( face.getID() ) ).at(1) );
             auto neighbor_cell_1_local_vertex_id_2 = static_cast< int64_t >( neighborCell1->getFaceLocalVertexToCellLocalVertexMaps().at( neighborCell1->getLocalFaceID( face.getID() ) ).at(2) );

             vertexdof::macroface::generated::apply_3D_macroface_vertexdof_to_vertexdof_replace(
                 dst_data,
                 src_data,
                 &src_data[offset_gl_0],
                 &src_data[offset_gl_1],
                 static_cast< int64_t >( level ),
                 neighbor_cell_0_local_vertex_id_0,
                 neighbor_cell_0_local_vertex_id_1,
                 neighbor_cell_0_local_vertex_id_2,
                 neighbor_cell_1_local_vertex_id_0,
                 neighbor_cell_1_local_vertex_id_1,
                 neighbor_cell_1_local_vertex_id_2,
                 opr_data );
           }
           else
           {
             vertexdof::macroface::apply3D< real_t >(
             level, face, *storage_, faceStencil3DID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
           }
        }
        else
        {
           if ( hhg::globalDefines::useGeneratedKernels )
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
  }

  this->timingTree_->stop( "Macro-Face" );

  this->timingTree_->start( "Macro-Cell" );

   for( const auto& it : storage_->getCells() )
   {
      Cell& cell = *it.second;

      const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
      if( testFlag( cellBC, flag ) )
      {
         if ( hhg::globalDefines::useGeneratedKernels )
         {
            auto opr_data = cell.getData( cellStencilID_ )->getData( level );
            real_t* src_data = cell.getData( src.getCellDataID() )->getPointer( level );
            real_t* dst_data = cell.getData( dst.getCellDataID() )->getPointer( level );
            if ( updateType == Replace )
            {
              vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_replace( dst_data, src_data, static_cast< int64_t >( level ), opr_data );
            }
            else if ( updateType == Add )
            {
              vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_add( dst_data, src_data, static_cast< int64_t >( level ), opr_data );
            }
         }
         else
         {
            vertexdof::macrocell::apply< real_t >( level, cell, cellStencilID_, src.getCellDataID(), dst.getCellDataID(), updateType );
         }

      }
   }

  this->timingTree_->stop( "Macro-Cell" );

   this->stopTiming( "Apply" );
}

template<class UFCOperator2D, class UFCOperator3D, bool Diagonal, bool Lumped, bool InvertDiagonal>
void P1ConstantOperator<UFCOperator2D, UFCOperator3D, Diagonal, Lumped, InvertDiagonal>::smooth_gs(
    const P1Function<real_t> &dst,
    const P1Function<real_t> &rhs,
    size_t level,
    DoFType flag) const
{
  smooth_sor( dst, rhs, 1.0, level, flag );
}


template < class UFCOperator2D, class UFCOperator3D, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< UFCOperator2D, UFCOperator3D, Diagonal, Lumped, InvertDiagonal >::smooth_sor(
    const P1Function< real_t >& dst,
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

   for( auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if( testFlag( vertexBC, flag ) )
      {
         vertexdof::macrovertex::smooth_sor(
                 vertex, vertexStencilID_, dst.getVertexDataID(), rhs.getVertexDataID(), level, relax );
      }
   }

   dst.communicate< Vertex, Edge >( level );

   for( auto& it : storage_->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if( testFlag( edgeBC, flag ) )
      {
         vertexdof::macroedge::smooth_sor< real_t >(
                 level, edge, edgeStencilID_, dst.getEdgeDataID(), rhs.getEdgeDataID(), relax );
      }
   }

   dst.communicate< Edge, Face >( level );

   for( auto& it : storage_->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if( testFlag( faceBC, flag ) )
      {
        if ( storage_->hasGlobalCells() )
        {
          vertexdof::macroface::smoothSOR3D< real_t >(
          level, face, *storage_, faceStencil3DID_, dst.getFaceDataID(), rhs.getFaceDataID(), relax );
        }
        else
        {
          vertexdof::macroface::smooth_sor< real_t >(
          level, face, faceStencilID_, dst.getFaceDataID(), rhs.getFaceDataID(), relax );
        }
      }
   }

   dst.communicate< Face, Cell >( level );

   for( auto& it : storage_->getCells() )
   {
      Cell& cell = *it.second;

      const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
      if( testFlag( cellBC, flag ) )
      {
         vertexdof::macrocell::smooth_sor< real_t >( level, cell, cellStencilID_, dst.getCellDataID(), rhs.getCellDataID(), relax );
      }
   }
}

template < class UFCOperator2D, class UFCOperator3D, bool Diagonal, bool Lumped, bool InvertDiagonal >
void P1ConstantOperator< UFCOperator2D, UFCOperator3D, Diagonal, Lumped, InvertDiagonal >::smooth_jac(
    const P1Function< real_t >& dst,
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

   for( auto& it : storage_->getVertices() )
   {
      Vertex& vertex = *it.second;

      const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
      if( testFlag( vertexBC, flag ) )
      {
         vertexdof::macrovertex::smooth_jac(
                 vertex, vertexStencilID_, dst.getVertexDataID(), rhs.getVertexDataID(), tmp.getVertexDataID(), level );
      }
   }

   for( auto& it : storage_->getEdges() )
   {
      Edge& edge = *it.second;

      const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
      if( testFlag( edgeBC, flag ) )
      {
         vertexdof::macroedge::smooth_jac< real_t >(
                 level, edge, edgeStencilID_, dst.getEdgeDataID(), rhs.getEdgeDataID(), tmp.getEdgeDataID() );
      }
   }

   for( auto& it : storage_->getFaces() )
   {
      Face& face = *it.second;

      const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
      if( testFlag( faceBC, flag ) )
      {
         vertexdof::macroface::smooth_jac< real_t >(
                 level, face, faceStencilID_, dst.getFaceDataID(), rhs.getFaceDataID(), tmp.getFaceDataID() );
      }
   }
}

template<class UFCOperator2D, class UFCOperator3D, bool Diagonal, bool Lumped, bool InvertDiagonal>
void P1ConstantOperator<UFCOperator2D, UFCOperator3D, Diagonal, Lumped, InvertDiagonal>::scale(real_t scalar) {
   for( uint_t level = minLevel_; level <= maxLevel_; ++level )
   {
      for( auto& it : storage_->getFaces() )
      {
         Face& face         = *it.second;
         auto  face_stencil = face.getData( faceStencilID_ )->getPointer( level );

         for( const auto& neighbor : vertexdof::macroface::neighborsWithCenter )
         {
            face_stencil[vertexdof::stencilIndexFromVertex( neighbor )] *= scalar;
         }
      }

      for( auto& it : storage_->getEdges() )
      {
         Edge& edge         = *it.second;
         auto  edge_stencil = edge.getData( edgeStencilID_ )->getPointer( level );

         edge_stencil[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )] *= scalar;

         for( const auto& neighbor : vertexdof::macroedge::neighborsOnEdgeFromVertexDoF )
         {
            edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] *= scalar;
         }

         for( const auto& neighbor : vertexdof::macroedge::neighborsOnSouthFaceFromVertexDoF )
         {
            edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] *= scalar;
         }

         if( edge.getNumNeighborFaces() == 2 )
         {
            for( const auto& neighbor : vertexdof::macroedge::neighborsOnNorthFaceFromVertexDoF )
            {
               edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] *= scalar;
            }
         }
      }

      for( auto& it : storage_->getVertices() )
      {
         Vertex& vertex         = *it.second;
         auto    vertex_stencil = vertex.getData( vertexStencilID_ )->getPointer( level );
         for( uint_t i = 0; i < vertex.getData( vertexStencilID_ )->getSize( level ); ++i )
         {
            vertex_stencil[i] *= scalar;
         }
      }
   }
}

template class P1ConstantOperator< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise >;

template class P1ConstantOperator< fenics::NoAssemble, fenics::NoAssemble >;
template class P1ConstantOperator< fenics::NoAssemble, fenics::UndefinedAssembly >;

template class P1ConstantOperator< p1_diffusion_cell_integral_0_otherwise, fenics::UndefinedAssembly, true >;

template class P1ConstantOperator< p1_stokes_epsilon_cell_integral_0_otherwise >;
template class P1ConstantOperator< p1_stokes_epsilon_cell_integral_1_otherwise >;
template class P1ConstantOperator< p1_stokes_epsilon_cell_integral_2_otherwise >;
template class P1ConstantOperator< p1_stokes_epsilon_cell_integral_3_otherwise >;

template class P1ConstantOperator< p1_div_cell_integral_0_otherwise, p1_tet_div_tet_cell_integral_0_otherwise >;
template class P1ConstantOperator< p1_div_cell_integral_1_otherwise, p1_tet_div_tet_cell_integral_1_otherwise >;
template class P1ConstantOperator< fenics::NoAssemble, p1_tet_div_tet_cell_integral_2_otherwise >;

template class P1ConstantOperator< p1_divt_cell_integral_0_otherwise, p1_tet_divt_tet_cell_integral_0_otherwise >;
template class P1ConstantOperator< p1_divt_cell_integral_1_otherwise, p1_tet_divt_tet_cell_integral_1_otherwise >;
template class P1ConstantOperator< fenics::NoAssemble, p1_tet_divt_tet_cell_integral_2_otherwise >;

template class P1ConstantOperator< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise >;
template class P1ConstantOperator< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise, false, true, true >;

template class P1ConstantOperator< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise >;

template class P1ConstantOperator< fenics::NoAssemble,        p2_tet_diffusion_cell_integral_0_otherwise >;
template class P1ConstantOperator< fenics::NoAssemble,        p2_tet_mass_cell_integral_0_otherwise >;
template class P1ConstantOperator< fenics::NoAssemble,        p2_tet_pspg_tet_cell_integral_0_otherwise >;

template class P1ConstantOperator< fenics::NoAssemble,        p2_tet_div_tet_cell_integral_0_otherwise >;
template class P1ConstantOperator< fenics::NoAssemble,        p2_tet_div_tet_cell_integral_1_otherwise >;
template class P1ConstantOperator< fenics::NoAssemble,        p2_tet_div_tet_cell_integral_2_otherwise >;

template class P1ConstantOperator< fenics::NoAssemble,        p2_tet_divt_tet_cell_integral_0_otherwise >;
template class P1ConstantOperator< fenics::NoAssemble,        p2_tet_divt_tet_cell_integral_1_otherwise >;
template class P1ConstantOperator< fenics::NoAssemble,        p2_tet_divt_tet_cell_integral_2_otherwise >;

template class P1ConstantOperator< fenics::NoAssemble,        p2_to_p1_tet_div_tet_cell_integral_0_otherwise >;
template class P1ConstantOperator< fenics::NoAssemble,        p2_to_p1_tet_div_tet_cell_integral_1_otherwise >;
template class P1ConstantOperator< fenics::UndefinedAssembly, p2_to_p1_tet_div_tet_cell_integral_2_otherwise >;
template class P1ConstantOperator< fenics::NoAssemble,        p2_to_p1_tet_div_tet_cell_integral_2_otherwise >;

template class P1ConstantOperator< fenics::NoAssemble,        p1_to_p2_tet_divt_tet_cell_integral_0_otherwise >;
template class P1ConstantOperator< fenics::NoAssemble,        p1_to_p2_tet_divt_tet_cell_integral_1_otherwise >;
template class P1ConstantOperator< fenics::UndefinedAssembly, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise >;
template class P1ConstantOperator< fenics::NoAssemble,        p1_to_p2_tet_divt_tet_cell_integral_2_otherwise >;

} // namespace hhg
