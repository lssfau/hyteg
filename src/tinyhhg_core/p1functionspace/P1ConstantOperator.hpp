#pragma once

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

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include <array>

#include "tinyhhg_core/Operator.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFIndexing.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/types/matrix.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

#include "P1DataHandling.hpp"
#include "P1Elements.hpp"
#include "P1Function.hpp"

namespace hhg {

using walberla::real_t;

template < class UFCOperator, bool Diagonal = false, bool Lumped = false, bool InvertDiagonal = false >
class P1ConstantOperator : public Operator< P1Function< real_t >, P1Function< real_t > >
{
 public:
   P1ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
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

      // Only assemble stencils if UFCOperator is specified
      if( !std::is_same< UFCOperator, fenics::NoAssemble >::value )
      {
         if ( storage_->hasGlobalCells() )
         {
            assembleStencils3D();
         }
         else
         {
            assembleStencils();
         }
      }
   }

   ~P1ConstantOperator() {}

   void scale( real_t scalar )
   {
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

   const PrimitiveDataID< StencilMemory< real_t >, Vertex >& getVertexStencilID() const { return vertexStencilID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Edge >& getEdgeStencilID() const { return edgeStencilID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Face >& getFaceStencilID() const { return faceStencilID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Cell >& getCellStencilID() const { return cellStencilID_; }

 private:
   void assembleStencils()
   {
      using namespace P1Elements::FaceVertexDoF;
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
               assembleP1LocalStencil(
                   P1GrayStencilMaps[i], P1GrayDoFMaps[i], local_stiffness_gray, face_stencil );
            }

            for( uint_t i = 0; i < P1BlueElements.size(); ++i )
            {
               assembleP1LocalStencil(
                   P1BlueStencilMaps[i], P1BlueDoFMaps[i], local_stiffness_blue, face_stencil );
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

            // first face
            Face* face = storage_->getFace( edge.neighborFaces()[0] );
            compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
            compute_local_stiffness( *face, level, local_stiffness_blue, fenics::BLUE );

            size_t start_id    = face->vertex_index( edge.neighborVertices()[0] );
            size_t end_id      = face->vertex_index( edge.neighborVertices()[1] );
            size_t opposite_id = face->vertex_index( face->get_vertex_opposite_to_edge( edge.getID() ) );

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
               // second face
               Face* face = storage_->getFace( edge.neighborFaces()[1] );
               compute_local_stiffness( *face, level, local_stiffness_gray, fenics::GRAY );
               compute_local_stiffness( *face, level, local_stiffness_blue, fenics::BLUE );

               size_t start_id    = face->vertex_index( edge.neighborVertices()[0] );
               size_t end_id      = face->vertex_index( edge.neighborVertices()[1] );
               size_t opposite_id = face->vertex_index( face->get_vertex_opposite_to_edge( edge.getID() ) );

               assembleP1LocalStencil( convertStencilDirectionsToIndices( elementNE ),
                                       {{start_id, end_id, opposite_id}},
                                       local_stiffness_gray,
                                       edge_stencil );
               assembleP1LocalStencil( convertStencilDirectionsToIndices( elementN ),
                                       {{opposite_id, start_id, end_id}},
                                       local_stiffness_blue,
                                       edge_stencil );
               assembleP1LocalStencil( convertStencilDirectionsToIndices( elementNW ),
                                       {{end_id, opposite_id, start_id}},
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

               std::array< uint_t, 3 > stencilMap;
               stencilMap[0] = 0;

               std::array< uint_t, 3 > dofMap;
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

   void assembleStencils3D()
   {
      for( uint_t level = minLevel_; level <= maxLevel_; level++ )
      {
         for ( const auto & it : storage_->getFaces() )
         {
            auto face = it.second;
            auto stencilSize   = face->getData( getFaceStencilID() )->getSize( level );
            auto stencilMemory = face->getData( getFaceStencilID() )->getPointer( level );
            UFCOperator ufcOperator;

            auto stencil = P1Elements::CellVertexDoF::assembleP1LocalStencil( storage_, *face, indexing::Index( 1, 1, 0 ), level, ufcOperator );

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
               const auto stencilIdx = vertexdof::stencilIndexFromVertex( stencilIt.first );
               stencilMemory[ stencilIdx ] = stencil[ stencilIt.first ];
            }
         }

         for ( const auto & it : storage_->getCells() )
         {
            auto cell          = it.second;
            auto stencilSize   = cell->getData( getCellStencilID() )->getSize( level );
            auto stencilMemory = cell->getData( getCellStencilID() )->getPointer( level );
            UFCOperator ufcOperator;

            auto stencil = P1Elements::CellVertexDoF::assembleP1LocalStencil( *cell, level, ufcOperator );

            for ( uint_t stencilEntryIdx = 0; stencilEntryIdx < stencilSize; stencilEntryIdx++ )
            {
               stencilMemory[ stencilEntryIdx ] = stencil[ stencilEntryIdx ];
            }
         }
      }
   }

 private:
   void apply_impl( P1Function< real_t >& src,
                    P1Function< real_t >& dst,
                    size_t                level,
                    DoFType               flag,
                    UpdateType            updateType = Replace )
   {
      src.getCommunicator( level )->communicate< Vertex, Edge >();
      src.getCommunicator( level )->communicate< Edge, Face >();
      src.getCommunicator( level )->communicate< Face, Cell >();

      // start pulling vertex halos
      src.getCommunicator( level )->startCommunication< Edge, Vertex >();

      // start pulling edge halos
      src.getCommunicator( level )->startCommunication< Face, Edge >();

      // start pulling face halos
      src.getCommunicator( level )->startCommunication< Cell, Face >();

      // end pulling vertex halos
      src.getCommunicator( level )->endCommunication< Edge, Vertex >();

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

      dst.getCommunicator( level )->startCommunication< Vertex, Edge >();

      // end pulling edge halos
      src.getCommunicator( level )->endCommunication< Face, Edge >();

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

      dst.getCommunicator( level )->endCommunication< Vertex, Edge >();

      dst.getCommunicator( level )->startCommunication< Edge, Face >();

      // end pulling face halos
      src.getCommunicator( level )->endCommunication< Cell, Face >();

      for( const auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if( testFlag( faceBC, flag ) )
         {
            vertexdof::macroface::apply< real_t >(
                level, face, faceStencilID_, src.getFaceDataID(), dst.getFaceDataID(), updateType );
         }
      }

      dst.getCommunicator( level )->endCommunication< Edge, Face >();

      for( const auto& it : storage_->getCells() )
      {
        Cell& cell = *it.second;

        const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );
        if ( testFlag( cellBC, flag ) )
        {
          vertexdof::macrocell::apply< real_t >(
          level, cell, cellStencilID_, src.getCellDataID(), dst.getCellDataID(), updateType );
        }
      }
   }

   void smooth_gs_impl( P1Function< real_t >& dst, P1Function< real_t >& rhs, size_t level, DoFType flag )
   {
      // start pulling vertex halos
      dst.getCommunicator( level )->startCommunication< Edge, Vertex >();

      // start pulling edge halos
      dst.getCommunicator( level )->startCommunication< Face, Edge >();

      // end pulling vertex halos
      dst.getCommunicator( level )->endCommunication< Edge, Vertex >();

      for( auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if( testFlag( vertexBC, flag ) )
         {
            vertexdof::macrovertex::smooth_gs< real_t >( vertex, vertexStencilID_, dst.getVertexDataID(), rhs.getVertexDataID(), level );
         }
      }

      dst.getCommunicator( level )->startCommunication< Vertex, Edge >();

      // end pulling edge halos
      dst.getCommunicator( level )->endCommunication< Face, Edge >();

      for( auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if( testFlag( edgeBC, flag ) )
         {
            vertexdof::macroedge::smooth_gs< real_t >( level, edge, edgeStencilID_, dst.getEdgeDataID(), rhs.getEdgeDataID() );
         }
      }

      dst.getCommunicator( level )->endCommunication< Vertex, Edge >();

      dst.getCommunicator( level )->startCommunication< Edge, Face >();

      for( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if( testFlag( faceBC, flag ) )
         {
            vertexdof::macroface::smooth_gs< real_t >( level, face, faceStencilID_, dst.getFaceDataID(), rhs.getFaceDataID() );
         }
      }

      dst.getCommunicator( level )->endCommunication< Edge, Face >();
   }

   void smooth_sor_impl( P1Function< real_t >& dst, P1Function< real_t >& rhs, real_t relax, size_t level, DoFType flag )
   {
      // start pulling vertex halos
      dst.getCommunicator( level )->startCommunication< Edge, Vertex >();

      // start pulling edge halos
      dst.getCommunicator( level )->startCommunication< Face, Edge >();

      // end pulling vertex halos
      dst.getCommunicator( level )->endCommunication< Edge, Vertex >();

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

      dst.getCommunicator( level )->startCommunication< Vertex, Edge >();

      // end pulling edge halos
      dst.getCommunicator( level )->endCommunication< Face, Edge >();

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

      dst.getCommunicator( level )->endCommunication< Vertex, Edge >();

      dst.getCommunicator( level )->startCommunication< Edge, Face >();

      for( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if( testFlag( faceBC, flag ) )
         {
            vertexdof::macroface::smooth_sor< real_t >(
                level, face, faceStencilID_, dst.getFaceDataID(), rhs.getFaceDataID(), relax );
         }
      }

      dst.getCommunicator( level )->endCommunication< Edge, Face >();
   }

   void smooth_jac_impl( P1Function< real_t >& dst,
                         P1Function< real_t >& rhs,
                         P1Function< real_t >& tmp,
                         size_t                level,
                         DoFType               flag )
   {
      // start pulling vertex halos
      tmp.getCommunicator( level )->startCommunication< Edge, Vertex >();

      // start pulling edge halos
      tmp.getCommunicator( level )->startCommunication< Face, Edge >();

      // end pulling vertex halos
      tmp.getCommunicator( level )->endCommunication< Edge, Vertex >();

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

      dst.getCommunicator( level )->startCommunication< Vertex, Edge >();

      // end pulling edge halos
      tmp.getCommunicator( level )->endCommunication< Face, Edge >();

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

      dst.getCommunicator( level )->endCommunication< Vertex, Edge >();

      dst.getCommunicator( level )->startCommunication< Edge, Face >();

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

      dst.getCommunicator( level )->endCommunication< Edge, Face >();
   }

   PrimitiveDataID< StencilMemory< real_t >, Vertex > vertexStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Edge >   edgeStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Face >   faceStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Cell >   cellStencilID_;

   void compute_local_stiffness( const Face& face, size_t level, Matrix3r& local_stiffness, fenics::ElementType element_type )
   {
      real_t coords[6];
      fenics::compute_micro_coords( face, level, coords, element_type );
      UFCOperator gen;
      gen.tabulate_tensor( local_stiffness.data(), NULL, coords, 0 );

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
};

typedef P1ConstantOperator< fenics::NoAssemble > P1ZeroOperator;

typedef P1ConstantOperator< p1_diffusion_cell_integral_0_otherwise >       P1LaplaceOperator;
typedef P1ConstantOperator< p1_diffusion_cell_integral_0_otherwise, true > P1DiagonalLaplaceOperator;

typedef P1ConstantOperator <p1_stokes_epsilon_cell_integral_0_otherwise > P1ConstantEpsilonOperator_11;
typedef P1ConstantOperator <p1_stokes_epsilon_cell_integral_1_otherwise > P1ConstantEpsilonOperator_12;
typedef P1ConstantOperator <p1_stokes_epsilon_cell_integral_2_otherwise > P1ConstantEpsilonOperator_21;
typedef P1ConstantOperator <p1_stokes_epsilon_cell_integral_3_otherwise > P1ConstantEpsilonOperator_22;

typedef P1ConstantOperator< p1_div_cell_integral_0_otherwise > P1DivxOperator;
typedef P1ConstantOperator< p1_div_cell_integral_1_otherwise > P1DivyOperator;

typedef P1ConstantOperator< p1_divt_cell_integral_0_otherwise > P1DivTxOperator;
typedef P1ConstantOperator< p1_divt_cell_integral_1_otherwise > P1DivTyOperator;

typedef P1ConstantOperator< p1_mass_cell_integral_0_otherwise >                    P1MassOperator;
typedef P1ConstantOperator< p1_mass_cell_integral_0_otherwise, false, true, true > P1LumpedInvMassOperator;

typedef P1ConstantOperator< p1_pspg_cell_integral_0_otherwise > P1PSPGOperator;

} // namespace hhg
