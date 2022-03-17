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
#pragma once



#include <array>

#include "core/OpenMP.h"

#include "hyteg/Stencil.hpp" 
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p1functionspace/variablestencil/VertexDoFVariableStencil.hpp"
#include "hyteg/solvers/Smoothables.hpp"

#include "P1Elements.hpp"

namespace hyteg {

using walberla::int_c;
using walberla::real_t;

template < class P1Form, bool Diagonal = false, bool Lumped = false, bool InvertDiagonal = false >
class P1Operator : public Operator< P1Function< real_t >, P1Function< real_t > >,
                   public GSSmoothable< P1Function< real_t > >,
                   public GSBackwardsSmoothable< P1Function< real_t > >,
                   public SORSmoothable< P1Function< real_t > >,
                   public SORBackwardsSmoothable< P1Function< real_t > >,
                   public WeightedJacobiSmoothable< P1Function< real_t > >,
                   public OperatorWithInverseDiagonal< P1Function< real_t > >
{
 protected:
   using Operator< P1Function< real_t >, P1Function< real_t > >::storage_;

 public:
   P1Operator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : P1Operator( storage, minLevel, maxLevel, P1Form() )
   {}

   P1Operator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel, const P1Form& form )
   : Operator( storage, minLevel, maxLevel )
   , form_( form )
   , formS_( form )
   , formN_( form )
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

      // pre-assemble vertex stencils
      if ( storage_->hasGlobalCells() )
      {
         assemble_stencil_vertices3D();
      }
      else
      {
         assemble_stencil_vertices();
      }
   }

   ~P1Operator() override = default;

   std::map< indexing::Index, vertexdof::macrocell::StencilMap_T > computeStencilsForCell( Cell& cell ) const
   {
      typedef stencilDirection sd;

      auto level = maxLevel_;

      auto& operatorData = cell.getData( cellStencilID_ )->getData( level );

      std::map< indexing::Index, vertexdof::macrocell::StencilMap_T > result;

      assemble_stencil_cell_init( cell, level );

      const uint_t rowsizeZ = levelinfo::num_microvertices_per_edge( level );
      uint_t       rowsizeY, rowsizeX;

      for ( uint_t k = 1; k < rowsizeZ - 3; ++k )
      {
         assemble_stencil_cell_init_z( k );

         rowsizeY = rowsizeZ - k;

         for ( uint_t j = 1; j < rowsizeY - 2; ++j )
         {
            assemble_stencil_cell_init_y( j );

            rowsizeX = rowsizeY - j;

            for ( uint_t i = 1; i < rowsizeX - 1; ++i )
            {
               indexing::Index idx{ i, j, k };

               assemble_stencil_cell( operatorData, i, j, k );

               result.insert_or_assign( idx, operatorData );
            }
         }
      }

      return result;
   }

   std::map< indexing::Index, std::vector< real_t > > computeStencilsForFace( Face& face ) const
   {
      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "Only available for 2D domains!" );
      }

      typedef stencilDirection sd;

      auto level = maxLevel_;

      // auto stencilSize = face.getData( faceStencilID_ )->getSize( level ); = 27
      uint_t stencilSize = 9;

      std::map< indexing::Index, std::vector< real_t > > result;

      uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );
      uint_t inner_rowsize = rowsize;

      real_t* opr_data = face.getData( faceStencilID_ )->getPointer( level );

      assemble_stencil_face_init( face, level );

      for ( uint_t j = 1; j < rowsize - 2; ++j )
      {
         assemble_stencil_face_init_y( j );

         for ( uint_t i = 1; i < inner_rowsize - 2; ++i )
         {
            indexing::Index idx{ i, j, 0 };

            assemble_stencil_face( opr_data, i, j );

            std::vector< real_t > stencil( stencilSize );

            std::copy( opr_data, opr_data + stencilSize, stencil.data() );

            result.insert_or_assign( idx, stencil );
         }

         --inner_rowsize;
      }

      return result;
   }

   void apply( const P1Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
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

         for ( int i = 0; i < int_c( edgeIDs.size() ); i++ )
         {
            Edge& edge = *this->getStorage()->getEdge( edgeIDs[uint_c( i )] );

            const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );

            if ( testFlag( edgeBC, flag ) )
            {
               apply_edge( edge, src.getEdgeDataID(), dst.getEdgeDataID(), level, updateType );
            }
         }
      }

      this->timingTree_->stop( "Macro-Edge" );

      this->timingTree_->start( "Macro-Face" );

      if ( level >= 2 )
      {
         std::vector< PrimitiveID > faceIDs = this->getStorage()->getFaceIDs();

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
                     apply_face3D_generated( face, src.getFaceDataID(), dst.getFaceDataID(), level, updateType );
                  }
                  else
                  {
                     apply_face3D( face, src.getFaceDataID(), dst.getFaceDataID(), level, updateType );
                  }
               }
               else
               {
                  if ( hyteg::globalDefines::useGeneratedKernels )
                  {
                     apply_face_generated( face, src.getFaceDataID(), dst.getFaceDataID(), level, updateType );
                  }
                  else
                  {
                     apply_face( face, src.getFaceDataID(), dst.getFaceDataID(), level, updateType );
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

         for ( int i = 0; i < int_c( cellIDs.size() ); i++ )
         {
            Cell& cell = *this->getStorage()->getCell( cellIDs[uint_c( i )] );

            const DoFType cellBC = dst.getBoundaryCondition().getBoundaryType( cell.getMeshBoundaryFlag() );

            if ( testFlag( cellBC, flag ) )
            {
               if ( hyteg::globalDefines::useGeneratedKernels )
               {
                  apply_cell_generated( cell, src.getCellDataID(), dst.getCellDataID(), level, updateType );
               }
               else
               {
                  apply_cell( cell, src.getCellDataID(), dst.getCellDataID(), level, updateType );
               }
            }
         }
      }

      this->timingTree_->stop( "Macro-Cell" );

      this->stopTiming( "Apply" );
   }

   void smooth_gs( const P1Function< real_t >& dst, const P1Function< real_t >& rhs, size_t level, DoFType flag ) const override

   {
      smooth_sor( dst, rhs, 1.0, level, flag );
   }

   void smooth_gs_backwards( const P1Function< real_t >& dst,
                             const P1Function< real_t >& rhs,
                             size_t                      level,
                             DoFType                     flag ) const override
   {
      smooth_sor_backwards( dst, rhs, 1.0, level, flag );
   }

   void smooth_sor( const P1Function< real_t >& dst,
                    const P1Function< real_t >& rhs,
                    real_t                      relax,
                    size_t                      level,
                    DoFType                     flag ) const override
   {
      smooth_sor( dst, rhs, relax, level, flag, false );
   }

   void smooth_sor( const P1Function< real_t >& dst,
                    const P1Function< real_t >& rhs,
                    real_t                      relax,
                    size_t                      level,
                    DoFType                     flag,
                    const bool&                 backwards ) const
   {
      if ( backwards )
      {
         if ( !backwards_sor_available() )
         {
            throw std::runtime_error( "Backward SOR not implemented for this operator." );
         }
         if ( !globalDefines::useGeneratedKernels )
         {
            throw std::runtime_error( "Backward SOR only implemented in generated kernels." );
         }

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

   void smooth_sor_backwards( const P1Function< real_t >& dst,
                              const P1Function< real_t >& rhs,
                              real_t                      relax,
                              size_t                      level,
                              DoFType                     flag ) const override
   {
      smooth_sor( dst, rhs, relax, level, flag, true );
   }

   void smooth_jac( const P1Function< real_t >& dst,
                    const P1Function< real_t >& rhs,
                    const P1Function< real_t >& src,
                    const real_t                relax,
                    size_t                      level,
                    DoFType                     flag ) const override
   {
      this->startTiming( "smooth_jac" );

      // compute the current residual
      this->apply( src, dst, level, flag );
      dst.assign( { real_c( 1 ), real_c( -1 ) }, { rhs, dst }, level, flag );

      // perform Jacobi update step
      dst.multElementwise( { *getInverseDiagonalValues(), dst }, level, flag );
      dst.assign( { 1.0, relax }, { src, dst }, level, flag );

      this->stopTiming( "smooth_jac" );
   }

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   // void computeDiagonalOperatorValues( bool use_variable_stencil_assembly = false )
   void computeDiagonalOperatorValues()
   {
      bool use_variable_stencil_assembly = false;
      computeDiagonalOperatorValues( false, use_variable_stencil_assembly );
   }

   /// Trigger (re)computation of inverse diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   // void computeInverseDiagonalOperatorValues( bool use_variable_stencil_assembly = false )
   void computeInverseDiagonalOperatorValues()
   {
      bool use_variable_stencil_assembly = false;
      computeDiagonalOperatorValues( true, use_variable_stencil_assembly );
   }

   std::shared_ptr< P1Function< real_t > > getDiagonalValues() const
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          diagonalValues_,
          "Diagonal values have not been assembled, call computeDiagonalOperatorValues() to set up this function." )
      return diagonalValues_;
   };

   std::shared_ptr< P1Function< real_t > > getInverseDiagonalValues() const override
   {
      WALBERLA_CHECK_NOT_NULLPTR(
          inverseDiagonalValues_,
          "Inverse diagonal values have not been assembled, call computeInverseDiagonalOperatorValues() to set up this function." )
      return inverseDiagonalValues_;
   };

   const PrimitiveDataID< StencilMemory< real_t >, Vertex >& getVertexStencilID() const { return vertexStencilID_; }

   const PrimitiveDataID< StencilMemory< real_t >, Edge >& getEdgeStencilID() const { return edgeStencilID_; }

   const PrimitiveDataID< LevelWiseMemory< vertexdof::macroface::StencilMap_T >, Edge >& getEdgeStencil3DID() const
   {
      return edgeStencil3DID_;
   }

   const PrimitiveDataID< StencilMemory< real_t >, Face >& getFaceStencilID() const { return faceStencilID_; }

   const PrimitiveDataID< LevelWiseMemory< vertexdof::macroface::StencilMap_T >, Face >& getFaceStencil3DID() const
   {
      return faceStencil3DID_;
   }

   const PrimitiveDataID< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell >& getCellStencilID() const
   {
      return cellStencilID_;
   }

 protected:
   // assemble stencils for all macro-vertices
   void assemble_stencil_vertices3D()
   {
      for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
      {
         for ( const auto& it : storage_->getVertices() )
         {
            auto vertex        = it.second;
            auto stencilSize   = vertex->getData( getVertexStencilID() )->getSize( level );
            auto stencilMemory = vertex->getData( getVertexStencilID() )->getPointer( level );

            form_.setGeometryMap( vertex->getGeometryMap() );
            auto stencil = P1Elements::P1Elements3D::assembleP1LocalStencil_new< P1Form >(
                storage_, *vertex, indexing::Index( 0, 0, 0 ), level, form_ );

            WALBERLA_ASSERT_EQUAL( stencilSize, stencil.size() );

            for ( uint_t i = 0; i < stencilSize; i++ )
            {
               stencilMemory[i] = stencil[i];
            }

            if constexpr ( Lumped )
            {
               for ( uint_t i = 1; i < stencilSize; i++ )
               {
                  stencilMemory[0] += stencilMemory[i];
                  stencilMemory[i] = 0;
               }
            }
            if constexpr ( Diagonal )
            {
               for ( uint_t i = 1; i < stencilSize; i++ )
               {
                  stencilMemory[i] = 0;
               }
            }
            if constexpr ( InvertDiagonal )
            {
               stencilMemory[0] = 1.0 / stencilMemory[0];
            }
         }
      }
   }

   // assemble stencils for all macro-vertices
   void assemble_stencil_vertices()
   {
      using namespace P1Elements::P1Elements2D;
      typedef stencilDirection sD;

      for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
      {
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
               d0 = ( face->coords[face->vertex_index(
                          storage_->getEdge( adj_edges[0] )->get_opposite_vertex( vertex.getID() ) )] -
                      x ) *
                    h;
               d2 = ( face->coords[face->vertex_index(
                          storage_->getEdge( adj_edges[1] )->get_opposite_vertex( vertex.getID() ) )] -
                      x ) *
                    h;

               Matrixr< 1, 3 > matrixRow;
               form_.integrateRow( 0, { { x, x + d0, x + d2 } }, matrixRow );

               uint_t i = 1;

               // iterate over adjacent edges
               for ( auto& edgeId : adj_edges )
               {
                  uint_t edge_idx = vertex.edge_index( edgeId ) + 1;

                  vertex_stencil[edge_idx] += matrixRow( 0, i );
                  i += 1;
               }

               // add contribution of center vertex
               vertex_stencil[0] += matrixRow( 0, 0 );

               ++neighborId;
            }

            if constexpr ( Lumped )
            {
               for ( uint_t i = 1; i < vertex.getData( vertexStencilID_ )->getSize( level ); ++i )
               {
                  vertex_stencil[0] += vertex_stencil[i];
                  vertex_stencil[i] = 0;
               }
            }

            if constexpr ( Diagonal )
            {
               for ( uint_t i = 1; i < vertex.getData( vertexStencilID_ )->getSize( level ); ++i )
               {
                  vertex_stencil[i] = 0;
               }
            }

            if constexpr ( InvertDiagonal )
            {
               vertex_stencil[0] = 1.0 / vertex_stencil[0];
            }
         }
      }
   }

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   ///
   /// \param invert if true, assembles the function carrying the inverse of the diagonal
   void computeDiagonalOperatorValues( bool invert, bool use_variable_stencil_assembly )
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
               auto stencilMemory = edge->getData( edgeStencilID_ )->getPointer( level );
               auto targetMemory  = edge->getData( targetFunction->getEdgeDataID() )->getPointer( level );

               if ( use_variable_stencil_assembly )
               {
                  assemble_variableStencil_edge_init( *edge, level );
               }
               else
               {
                  assemble_stencil_edge_init( *edge, level );
               }

               for ( auto idx : vertexdof::macroedge::Iterator( level ) )
               {
                  if ( variableStencil() )
                  {
                     if ( use_variable_stencil_assembly )
                     {
                        assemble_variableStencil_edge( stencilMemory, idx.x() );
                     }
                     else
                     {
                        assemble_stencil_edge( stencilMemory, idx.x() );
                     }
                  }

                  targetMemory[vertexdof::macroedge::index( level, idx.x() )] =
                      stencilMemory[vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C )];
               }
            }
         }

         if ( level >= 1 )
         {
            for ( const auto& it : storage_->getFaces() )
            {
               auto  face          = it.second;
               auto  stencilMemory = face->getData( faceStencilID_ )->getPointer( level );
               auto& stencilMap    = face->getData( faceStencil3DID_ )->getData( level );
               auto  targetMemory  = face->getData( targetFunction->getFaceDataID() )->getPointer( level );

               if ( use_variable_stencil_assembly )
               {
                  assemble_variableStencil_face_init( *face, level );
               }
               else
               {
                  assemble_stencil_face_init( *face, level );
               }

               const uint_t rowsizeY = levelinfo::num_microvertices_per_edge( level );
               uint_t       rowsizeX;

               for ( uint_t j = 1; j < rowsizeY - 2; ++j )
               {
                  if ( !use_variable_stencil_assembly )
                  {
                     assemble_stencil_face_init_y( j );
                  }

                  rowsizeX = rowsizeY - j;

                  for ( uint_t i = 1; i < rowsizeX - 1; ++i )
                  {
                     real_t centerValue = 0;

                     if ( storage_->hasGlobalCells() )
                     {
                        if ( variableStencil() )
                        {
                           if ( use_variable_stencil_assembly )
                           {
                              assemble_variableStencil_face3D( stencilMap, i, j );
                           }
                           else
                           {
                              assemble_stencil_face3D( stencilMap, i, j );
                           }
                        }

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
                        if ( variableStencil() )
                        {
                           if ( use_variable_stencil_assembly )
                           {
                              assemble_variableStencil_face( stencilMemory, i, j );
                           }
                           else
                           {
                              assemble_stencil_face( stencilMemory, i, j );
                           }
                        }

                        centerValue = stencilMemory[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];
                     }

                     targetMemory[vertexdof::macroface::index( level, i, j )] = centerValue;
                  }
               }

               // for (auto idx : vertexdof::macroface::Iterator(level))
               // {
               //    real_t centerValue = 0;

               //    if (storage_->hasGlobalCells())
               //    {
               //       if (variableStencil())
               //       {
               //          assemble_stencil_face3D(stencilMap, idx.x(), idx.y());
               //       }

               //       for (uint_t neighborCellID = 0; neighborCellID < face->getNumNeighborCells(); neighborCellID++)
               //       {
               //          for (auto stencilIt : stencilMap[neighborCellID])
               //          {
               //             if (stencilIt.first == indexing::IndexIncrement({0, 0, 0}))
               //             {
               //                centerValue += stencilIt.second;
               //             }
               //          }
               //       }
               //    }
               //    else
               //    {
               //       if (variableStencil())
               //       {
               //          face->getData(faceStencilID_)->setToZero(level);
               //          assemble_stencil_face(stencilMemory, idx.x(), idx.y());
               //       }

               //       centerValue = stencilMemory[vertexdof::stencilIndexFromVertex(stencilDirection::VERTEX_C)];
               //    }

               //    targetMemory[vertexdof::macroface::index(level, idx.x(), idx.y())] = centerValue;
               // }
            }
         }

         if ( level >= 2 )
         {
            for ( const auto& it : storage_->getCells() )
            {
               auto  cell         = it.second;
               auto& stencilMap   = cell->getData( cellStencilID_ )->getData( level );
               auto  targetMemory = cell->getData( targetFunction->getCellDataID() )->getPointer( level );

               if ( use_variable_stencil_assembly )
               {
                  assemble_variableStencil_cell_init( *cell, level );
               }
               else
               {
                  assemble_stencil_cell_init( *cell, level );
               }

               const uint_t rowsizeZ = levelinfo::num_microvertices_per_edge( level );
               uint_t       rowsizeY, rowsizeX;

               for ( uint_t k = 1; k < rowsizeZ - 3; ++k )
               {
                  if ( !use_variable_stencil_assembly )
                  {
                     assemble_stencil_cell_init_z( k );
                  }

                  rowsizeY = rowsizeZ - k;

                  for ( uint_t j = 1; j < rowsizeY - 2; ++j )
                  {
                     if ( !use_variable_stencil_assembly )
                     {
                        assemble_stencil_cell_init_y( j );
                     }

                     rowsizeX = rowsizeY - j;

                     for ( uint_t i = 1; i < rowsizeX - 1; ++i )
                     {
                        if ( variableStencil() )
                        {
                           if ( use_variable_stencil_assembly )
                           {
                              assemble_variableStencil_cell( stencilMap, i, j, k );
                           }
                           else
                           {
                              assemble_stencil_cell( stencilMap, i, j, k );
                           }
                        }

                        real_t centerValue = stencilMap[indexing::IndexIncrement( { 0, 0, 0 } )];

                        targetMemory[vertexdof::macrocell::index( level, i, j, k )] = centerValue;
                     }
                  }
               }

               // for (auto idx : vertexdof::macrocell::Iterator(level))
               // {
               //    if (variableStencil())
               //    {
               //       assemble_stencil_cell(stencilMap, idx.x(), idx.y(), idx.z());
               //    }

               //    real_t centerValue = stencilMap[indexing::IndexIncrement({0, 0, 0})];

               //    targetMemory[vertexdof::macrocell::index(level, idx.x(), idx.y(), idx.z())] = centerValue;
               // }
            }
         }

         if ( invert )
         {
            targetFunction->invertElementwise( level, All, false );
         }
      }
   }

   void smooth_sor_macro_vertices( const P1Function< real_t >& dst,
                                   const P1Function< real_t >& rhs,
                                   real_t                      relax,
                                   size_t                      level,
                                   DoFType                     flag,
                                   const bool&                 backwards = false ) const
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

   void smooth_sor_macro_edges( const P1Function< real_t >& dst,
                                const P1Function< real_t >& rhs,
                                real_t                      relax,
                                size_t                      level,
                                DoFType                     flag,
                                const bool&                 backwards = false ) const
   {
      this->timingTree_->start( "Macro-Edge" );

      for ( auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );

         if ( testFlag( edgeBC, flag ) )
         {
            smooth_sor_edge( edge, dst.getEdgeDataID(), rhs.getEdgeDataID(), level, relax, backwards );
         }
      }

      this->timingTree_->stop( "Macro-Edge" );
   }

   void smooth_sor_macro_faces( const P1Function< real_t >& dst,
                                const P1Function< real_t >& rhs,
                                real_t                      relax,
                                size_t                      level,
                                DoFType                     flag,
                                const bool&                 backwards = false ) const
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
                  smooth_sor_face3D_generated( face, dst.getFaceDataID(), rhs.getFaceDataID(), level, relax, backwards );
               }
               else
               {
                  smooth_sor_face3D( face, dst.getFaceDataID(), rhs.getFaceDataID(), level, relax, backwards );
               }
            }
            else
            {
               if ( globalDefines::useGeneratedKernels )
               {
                  smooth_sor_face_generated( face, dst.getFaceDataID(), rhs.getFaceDataID(), level, relax, backwards );
               }
               else
               {
                  smooth_sor_face( face, dst.getFaceDataID(), rhs.getFaceDataID(), level, relax, backwards );
               }
            }
         }
      }

      this->timingTree_->stop( "Macro-Face" );
   }

   void smooth_sor_macro_cells( const P1Function< real_t >& dst,
                                const P1Function< real_t >& rhs,
                                real_t                      relax,
                                size_t                      level,
                                DoFType                     flag,
                                const bool&                 backwards = false ) const
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
               smooth_sor_cell_generated( cell, dst.getCellDataID(), rhs.getCellDataID(), level, relax, backwards );
            }
            else
            {
               smooth_sor_cell( cell, dst.getCellDataID(), rhs.getCellDataID(), level, relax, backwards );
            }
         }
      }

      this->timingTree_->stop( "Macro-Cell" );
   }

   // apply the operator to all DoF on a given macro-edge
   inline void apply_edge( Edge&                                                    edge,
                           const PrimitiveDataID< FunctionMemory< real_t >, Edge >& srcId,
                           const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstId,
                           const uint_t&                                            level,
                           UpdateType                                               update ) const
   {
      using sD       = stencilDirection;
      size_t rowsize = levelinfo::num_microvertices_per_edge( level );

      auto opr_data = edge.getData( edgeStencilID_ )->getPointer( level );
      auto src      = edge.getData( srcId )->getPointer( level );
      auto dst      = edge.getData( dstId )->getPointer( level );

      assemble_stencil_edge_init( edge, level );

      real_t tmp;

      for ( size_t i = 1; i < rowsize - 1; ++i )
      {
         if ( variableStencil() )
         {
            assemble_stencil_edge( opr_data, i );
         }

         const auto stencilIdxW = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_W );
         const auto stencilIdxC = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_C );
         const auto stencilIdxE = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_E );

         const auto dofIdxW = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_W );
         const auto dofIdxC = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C );
         const auto dofIdxE = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_E );

         tmp = opr_data[stencilIdxW] * src[dofIdxW] + opr_data[stencilIdxC] * src[dofIdxC] + opr_data[stencilIdxE] * src[dofIdxE];

         for ( uint_t neighborFace = 0; neighborFace < edge.getNumNeighborFaces(); neighborFace++ )
         {
            const auto stencilIdxWNeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_W, neighborFace );
            const auto stencilIdxENeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_E, neighborFace );
            const auto stencilWeightW          = opr_data[stencilIdxWNeighborFace];
            const auto stencilWeightE          = opr_data[stencilIdxENeighborFace];
            const auto dofIdxWNeighborFace =
                vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_W );
            const auto dofIdxENeighborFace =
                vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_E );
            tmp += stencilWeightW * src[dofIdxWNeighborFace] + stencilWeightE * src[dofIdxENeighborFace];
         }

         for ( uint_t neighborCell = 0; neighborCell < edge.getNumNeighborCells(); neighborCell++ )
         {
            const auto stencilIdx = vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell, edge.getNumNeighborFaces() );
            const auto dofIdx =
                vertexdof::macroedge::indexFromVertexOnNeighborCell( level, i, neighborCell, edge.getNumNeighborFaces() );
            tmp += opr_data[stencilIdx] * src[dofIdx];
         }

         if ( update == Replace )
         {
            dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] = tmp;
         }
         else if ( update == Add )
         {
            dst[vertexdof::macroedge::indexFromVertex( level, i, stencilDirection::VERTEX_C )] += tmp;
         }
      }
   }

   // apply the operator to all DoF on a given macro-face
   inline void apply_face3D( Face&                                                    face,
                             const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcId,
                             const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                             const uint_t&                                            level,
                             UpdateType                                               update ) const
   {
      auto&   opr_data = face.getData( faceStencil3DID_ )->getData( level );
      real_t* src      = face.getData( srcId )->getPointer( level );
      real_t* dst      = face.getData( dstId )->getPointer( level );

      assemble_stencil_face_init( face, level );

      // todo loop ij
      for ( const auto& idxIt : vertexdof::macroface::Iterator( level, 1 ) )
      {
         if ( variableStencil() )
         {
            assemble_stencil_face3D( opr_data, idxIt.x(), idxIt.y() );
         }

         real_t tmp = real_c( 0 );

         for ( uint_t neighborCellIdx = 0; neighborCellIdx < face.getNumNeighborCells(); neighborCellIdx++ )
         {
            auto neighborCell = storage_->getCell( face.neighborCells().at( neighborCellIdx ) );
            auto centerIndexInCell =
                vertexdof::macroface::getIndexInNeighboringMacroCell( idxIt, face, neighborCellIdx, *storage_, level );

            for ( auto stencilIt : opr_data[neighborCellIdx] )
            {
               auto weight               = stencilIt.second;
               auto leafIndexInMacroCell = centerIndexInCell + stencilIt.first;
               auto leafIndexInMacroFace = vertexdof::macrocell::getIndexInNeighboringMacroFace(
                   leafIndexInMacroCell, *neighborCell, neighborCell->getLocalFaceID( face.getID() ), *storage_, level );

               uint_t leafArrayIndexInMacroFace;

               if ( leafIndexInMacroFace.z() == 0 )
               {
                  leafArrayIndexInMacroFace =
                      vertexdof::macroface::index( level, leafIndexInMacroFace.x(), leafIndexInMacroFace.y() );
               }
               else
               {
                  WALBERLA_ASSERT_EQUAL( leafIndexInMacroFace.z(), 1 );
                  leafArrayIndexInMacroFace =
                      vertexdof::macroface::index( level, leafIndexInMacroFace.x(), leafIndexInMacroFace.y(), neighborCellIdx );
               }

               tmp += weight * src[leafArrayIndexInMacroFace];
            }
         }

         if ( update == Replace )
         {
            dst[vertexdof::macroface::index( level, idxIt.x(), idxIt.y() )] = tmp;
         }
         else if ( update == Add )
         {
            dst[vertexdof::macroface::index( level, idxIt.x(), idxIt.y() )] += tmp;
         }
      }
   }

   // apply the operator to all DoF on a given macro-face
   inline void apply_face( Face&                                                    face,
                           const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcId,
                           const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                           const uint_t&                                            level,
                           UpdateType                                               update ) const
   {
      uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );
      uint_t inner_rowsize = rowsize;

      real_t* opr_data = face.getData( faceStencilID_ )->getPointer( level );
      real_t* src      = face.getData( srcId )->getPointer( level );
      real_t* dst      = face.getData( dstId )->getPointer( level );

      assemble_stencil_face_init( face, level );

      real_t tmp = real_c( 0 );

      for ( uint_t j = 1; j < rowsize - 2; ++j )
      {
         assemble_stencil_face_init_y( j );

         for ( uint_t i = 1; i < inner_rowsize - 2; ++i )
         {
            if ( variableStencil() )
            {
               assemble_stencil_face( opr_data, i, j );
            }

            if ( face.getNumNeighborCells() == 0 )
            {
               static_assert( vertexdof::macroface::neighborsWithoutCenter.size() == 6, "Neighbors array has wrong size" );
               tmp = real_c( 0 );

               for ( const auto direction : vertexdof::macroface::neighborsWithCenter )
               {
                  tmp += opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                         src[vertexdof::macroface::indexFromVertex( level, i, j, direction )];
               }
            }
            else if ( face.getNumNeighborCells() == 1 )
            {
               tmp = real_c( 0 );

               for ( const auto direction : vertexdof::macroface::neighborsWithOneNeighborCellWithCenter )
               {
                  tmp += opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                         src[vertexdof::macroface::indexFromVertex( level, i, j, direction )];
               }
            }
            else if ( face.getNumNeighborCells() == 2 )
            {
               tmp = real_c( 0 );

               for ( const auto direction : vertexdof::macroface::neighborsWithTwoNeighborCellsWithCenter )
               {
                  tmp += opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                         src[vertexdof::macroface::indexFromVertex( level, i, j, direction )];
               }
            }

            WALBERLA_ASSERT_LESS( face.getNumNeighborCells(), 3 );

            if ( update == Replace )
            {
               dst[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] = tmp;
            }
            else
            {
               dst[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] += tmp;
            }
         }

         --inner_rowsize;
      }
   }

   // apply the operator to all DoF on a given macro-cell
   inline void apply_cell( Cell&                                                    cell,
                           const PrimitiveDataID< FunctionMemory< real_t >, Cell >& srcId,
                           const PrimitiveDataID< FunctionMemory< real_t >, Cell >& dstId,
                           const uint_t&                                            level,
                           UpdateType                                               update ) const
   {
      typedef stencilDirection sd;

      auto&         operatorData = cell.getData( cellStencilID_ )->getData( level );
      const real_t* src          = cell.getData( srcId )->getPointer( level );
      real_t*       dst          = cell.getData( dstId )->getPointer( level );

      assemble_stencil_cell_init( cell, level );

      real_t tmp;

      const uint_t rowsizeZ = levelinfo::num_microvertices_per_edge( level );
      uint_t       rowsizeY, rowsizeX;

      // skip level 0 (no interior points)
      if ( level == 0 )
         return;

      for ( uint_t k = 1; k < rowsizeZ - 3; ++k )
      {
         assemble_stencil_cell_init_z( k );

         rowsizeY = rowsizeZ - k;

         for ( uint_t j = 1; j < rowsizeY - 2; ++j )
         {
            assemble_stencil_cell_init_y( j );

            rowsizeX = rowsizeY - j;

            for ( uint_t i = 1; i < rowsizeX - 1; ++i )
            {
               if ( variableStencil() )
               {
                  assemble_stencil_cell( operatorData, i, j, k );
               }

               const uint_t centerIdx = vertexdof::macrocell::indexFromVertex( level, i, j, k, sd::VERTEX_C );

               tmp = operatorData.at( { 0, 0, 0 } ) * src[centerIdx];

               for ( const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter )
               {
                  const uint_t idx = vertexdof::macrocell::indexFromVertex( level, i, j, k, neighbor );
                  tmp += operatorData[vertexdof::logicalIndexOffsetFromVertex( neighbor )] * src[idx];
               }

               if ( update == Replace )
               {
                  dst[centerIdx] = tmp;
               }
               else
               {
                  dst[centerIdx] += tmp;
               }
            }
         }
      }

      // for (const auto& it : vertexdof::macrocell::Iterator(level, 1))
      // {
      //    const uint_t x = it.x();
      //    const uint_t y = it.y();
      //    const uint_t z = it.z();

      //    if (variableStencil())
      //    {
      //       assemble_stencil_cell(operatorData, x, y, z);
      //    }

      //    const uint_t centerIdx = vertexdof::macrocell::indexFromVertex(level, x, y, z, sd::VERTEX_C);

      //    tmp = operatorData.at({0, 0, 0}) * src[ centerIdx ];

      //    for (const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter)
      //    {
      //       const uint_t idx        = vertexdof::macrocell::indexFromVertex(level, x, y, z, neighbor);
      //       WALBERLA_ASSERT_GREATER(operatorData.count(vertexdof::logicalIndexOffsetFromVertex(neighbor)), 0);
      //       tmp += operatorData.at(vertexdof::logicalIndexOffsetFromVertex(neighbor)) * src[ idx ];
      //    }

      //    if (update == Replace)
      //    {
      //       dst[ centerIdx ] = tmp;
      //    }
      //    else
      //    {
      //       dst[ centerIdx ] += tmp;

      //    }
      // }
   }

   // apply the sor-operator to all DoF on a given macro-edge
   inline void smooth_sor_edge( Edge&                                                    edge,
                                const PrimitiveDataID< FunctionMemory< real_t >, Edge >& dstId,
                                const PrimitiveDataID< FunctionMemory< real_t >, Edge >& rhsId,
                                const uint_t&                                            level,
                                real_t                                                   relax,
                                const bool&                                              backwards = false ) const
   {
      using sD       = stencilDirection;
      size_t rowsize = levelinfo::num_microvertices_per_edge( level );

      auto opr_data = edge.getData( edgeStencilID_ )->getPointer( level );
      auto rhs      = edge.getData( rhsId )->getPointer( level );
      auto dst      = edge.getData( dstId )->getPointer( level );

      const auto stencilIdxW = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_W );
      const auto stencilIdxC = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_C );
      const auto stencilIdxE = vertexdof::macroedge::stencilIndexOnEdge( sD::VERTEX_E );

      assemble_stencil_edge_init( edge, level );

      auto invCenterWeight = 1.0 / opr_data[stencilIdxC];

      real_t tmp;

      const int start = backwards ? (int) rowsize - 2 : 1;
      const int stop  = backwards ? 0 : (int) rowsize - 1;
      const int incr  = backwards ? -1 : 1;

      for ( int ii = start; ii != stop; ii += incr )
      {
         const uint_t i = uint_c( ii );

         if ( variableStencil() )
         {
            assemble_stencil_edge( opr_data, i );
            invCenterWeight = 1.0 / opr_data[stencilIdxC];
         }

         const auto dofIdxW = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_W );
         const auto dofIdxC = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_C );
         const auto dofIdxE = vertexdof::macroedge::indexFromVertex( level, i, sD::VERTEX_E );

         tmp = rhs[dofIdxC];

         tmp -= opr_data[stencilIdxW] * dst[dofIdxW] + opr_data[stencilIdxE] * dst[dofIdxE];

         for ( uint_t neighborFace = 0; neighborFace < edge.getNumNeighborFaces(); neighborFace++ )
         {
            const auto stencilIdxWNeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_W, neighborFace );
            const auto stencilIdxENeighborFace = vertexdof::macroedge::stencilIndexOnNeighborFace( sD::VERTEX_E, neighborFace );
            const auto stencilWeightW          = opr_data[stencilIdxWNeighborFace];
            const auto stencilWeightE          = opr_data[stencilIdxENeighborFace];
            const auto dofIdxWNeighborFace =
                vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_W );
            const auto dofIdxENeighborFace =
                vertexdof::macroedge::indexFromVertexOnNeighborFace( level, i, neighborFace, sD::VERTEX_E );
            tmp -= stencilWeightW * dst[dofIdxWNeighborFace] + stencilWeightE * dst[dofIdxENeighborFace];
         }

         for ( uint_t neighborCell = 0; neighborCell < edge.getNumNeighborCells(); neighborCell++ )
         {
            const auto stencilIdx = vertexdof::macroedge::stencilIndexOnNeighborCell( neighborCell, edge.getNumNeighborFaces() );
            const auto dofIdx =
                vertexdof::macroedge::indexFromVertexOnNeighborCell( level, i, neighborCell, edge.getNumNeighborFaces() );
            tmp -= opr_data[stencilIdx] * dst[dofIdx];
         }

         dst[dofIdxC] = ( 1.0 - relax ) * dst[dofIdxC] + relax * invCenterWeight * tmp;
      }
   }

   // apply the sor-operator to all DoF on a given macro-face
   inline void smooth_sor_face3D( Face&                                                    face,
                                  const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                                  const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsId,
                                  const uint_t&                                            level,
                                  real_t                                                   relax,
                                  const bool&                                              backwards = false ) const
   {
      WALBERLA_UNUSED( backwards );

      auto& opr_data = face.getData( faceStencil3DID_ )->getData( level );
      auto  rhs      = face.getData( rhsId )->getPointer( level );
      auto  dst      = face.getData( dstId )->getPointer( level );

      assemble_stencil_face_init( face, level );

      real_t centerWeight = real_c( 0 );

      for ( uint_t neighborCellIdx = 0; neighborCellIdx < face.getNumNeighborCells(); neighborCellIdx++ )
      {
         centerWeight += opr_data[neighborCellIdx][{ 0, 0, 0 }];
      }

      auto invCenterWeight = 1.0 / centerWeight;

      // todo loop ij
      for ( const auto& idxIt : vertexdof::macroface::Iterator( level, 1 ) )
      {
         real_t tmp = rhs[vertexdof::macroface::index( level, idxIt.x(), idxIt.y() )];

         if ( variableStencil() )
         {
            assemble_stencil_face3D( opr_data, idxIt.x(), idxIt.y() );
            centerWeight = real_c( 0 );

            for ( uint_t neighborCellIdx = 0; neighborCellIdx < face.getNumNeighborCells(); neighborCellIdx++ )
            {
               centerWeight += opr_data[neighborCellIdx][{ 0, 0, 0 }];
            }

            invCenterWeight = 1.0 / centerWeight;
         }

         for ( uint_t neighborCellIdx = 0; neighborCellIdx < face.getNumNeighborCells(); neighborCellIdx++ )
         {
            auto neighborCell = storage_->getCell( face.neighborCells().at( neighborCellIdx ) );
            auto centerIndexInCell =
                vertexdof::macroface::getIndexInNeighboringMacroCell( idxIt, face, neighborCellIdx, *storage_, level );

            for ( auto stencilIt : opr_data[neighborCellIdx] )
            {
               if ( stencilIt.first == indexing::IndexIncrement( { 0, 0, 0 } ) )
                  continue;

               auto weight               = stencilIt.second;
               auto leafIndexInMacroCell = centerIndexInCell + stencilIt.first;
               auto leafIndexInMacroFace = vertexdof::macrocell::getIndexInNeighboringMacroFace(
                   leafIndexInMacroCell, *neighborCell, neighborCell->getLocalFaceID( face.getID() ), *storage_, level );

               uint_t leafArrayIndexInMacroFace;

               if ( leafIndexInMacroFace.z() == 0 )
               {
                  leafArrayIndexInMacroFace =
                      vertexdof::macroface::index( level, leafIndexInMacroFace.x(), leafIndexInMacroFace.y() );
               }
               else
               {
                  WALBERLA_ASSERT_EQUAL( leafIndexInMacroFace.z(), 1 );
                  leafArrayIndexInMacroFace =
                      vertexdof::macroface::index( level, leafIndexInMacroFace.x(), leafIndexInMacroFace.y(), neighborCellIdx );
               }

               tmp -= weight * dst[leafArrayIndexInMacroFace];
            }
         }

         dst[vertexdof::macroface::index( level, idxIt.x(), idxIt.y() )] =
             ( 1.0 - relax ) * dst[vertexdof::macroface::index( level, idxIt.x(), idxIt.y() )] + relax * tmp * invCenterWeight;
      }
   }

   // apply the sor-operator to all DoF on a given macro-face
   inline void smooth_sor_face( Face&                                                    face,
                                const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                                const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsId,
                                const uint_t&                                            level,
                                real_t                                                   relax,
                                const bool&                                              backwards = false ) const
   {
      WALBERLA_UNUSED( backwards );

      uint_t rowsize       = levelinfo::num_microvertices_per_edge( level );
      uint_t inner_rowsize = rowsize;

      auto opr_data = face.getData( faceStencilID_ )->getPointer( level );
      auto dst      = face.getData( dstId )->getPointer( level );
      auto rhs      = face.getData( rhsId )->getPointer( level );

      assemble_stencil_face_init( face, level );

      auto invCenterWeight = 1.0 / opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];

      real_t tmp;

      for ( uint_t j = 1; j < rowsize - 2; ++j )
      {
         assemble_stencil_face_init_y( j );

         for ( uint_t i = 1; i < inner_rowsize - 2; ++i )
         {
            if ( variableStencil() )
            {
               assemble_stencil_face( opr_data, i, j );
               invCenterWeight = 1.0 / opr_data[vertexdof::stencilIndexFromVertex( stencilDirection::VERTEX_C )];
            }

            tmp = rhs[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )];

            if ( face.getNumNeighborCells() == 0 )
            {
               for ( const auto direction : vertexdof::macroface::neighborsWithoutCenter )
               {
                  tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                         dst[vertexdof::macroface::indexFromVertex( level, i, j, direction )];
               }
            }
            else if ( face.getNumNeighborCells() == 1 )
            {
               for ( const auto direction : vertexdof::macroface::neighborsWithOneNeighborCellWithoutCenter )
               {
                  tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                         dst[vertexdof::macroface::indexFromVertex( level, i, j, direction )];
               }
            }
            else if ( face.getNumNeighborCells() == 2 )
            {
               for ( const auto direction : vertexdof::macroface::neighborsWithTwoNeighborCellsWithoutCenter )
               {
                  tmp -= opr_data[vertexdof::stencilIndexFromVertex( direction )] *
                         dst[vertexdof::macroface::indexFromVertex( level, i, j, direction )];
               }
            }

            WALBERLA_ASSERT_LESS( face.getNumNeighborCells(), 3 );

            dst[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] =
                ( 1.0 - relax ) * dst[vertexdof::macroface::indexFromVertex( level, i, j, stencilDirection::VERTEX_C )] +
                relax * tmp * invCenterWeight;
         }

         --inner_rowsize;
      }
   }

   // apply the sor-operator to all DoF on a given macro-cell
   virtual void smooth_sor_cell( Cell&                                                    cell,
                                 const PrimitiveDataID< FunctionMemory< real_t >, Cell >& dstId,
                                 const PrimitiveDataID< FunctionMemory< real_t >, Cell >& rhsId,
                                 const uint_t&                                            level,
                                 real_t                                                   relax,
                                 const bool&                                              backwards = false ) const
   {
      WALBERLA_UNUSED( backwards );

      typedef stencilDirection sd;

      auto&         operatorData = cell.getData( cellStencilID_ )->getData( level );
      const real_t* rhs          = cell.getData( rhsId )->getPointer( level );
      real_t*       dst          = cell.getData( dstId )->getPointer( level );

      assemble_stencil_cell_init( cell, level );

      real_t tmp;

      auto inverseCenterWeight = 1.0 / operatorData[{ 0, 0, 0 }];

      const uint_t rowsizeZ = levelinfo::num_microvertices_per_edge( level );
      uint_t       rowsizeY, rowsizeX;

      // skip level 0 (no interior points)
      if ( level == 0 )
         return;

      for ( uint_t k = 1; k < rowsizeZ - 3; ++k )
      {
         assemble_stencil_cell_init_z( k );

         rowsizeY = rowsizeZ - k;

         for ( uint_t j = 1; j < rowsizeY - 2; ++j )
         {
            assemble_stencil_cell_init_y( j );

            rowsizeX = rowsizeY - j;

            for ( uint_t i = 1; i < rowsizeX - 1; ++i )
            {
               if ( variableStencil() )
               {
                  assemble_stencil_cell( operatorData, i, j, k );
                  inverseCenterWeight = 1.0 / operatorData[{ 0, 0, 0 }];
               }

               const uint_t centerIdx = vertexdof::macrocell::indexFromVertex( level, i, j, k, sd::VERTEX_C );

               tmp = rhs[centerIdx];

               for ( const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter )
               {
                  const uint_t idx = vertexdof::macrocell::indexFromVertex( level, i, j, k, neighbor );
                  tmp -= operatorData[vertexdof::logicalIndexOffsetFromVertex( neighbor )] * dst[idx];
               }

               dst[centerIdx] = ( 1.0 - relax ) * dst[centerIdx] + tmp * relax * inverseCenterWeight;
            }
         }
      }

      // for (const auto& it : vertexdof::macrocell::Iterator(level, 1))
      // {
      //    const uint_t x = it.x();
      //    const uint_t y = it.y();
      //    const uint_t z = it.z();

      //    if (variableStencil())
      //    {
      //       assemble_stencil_cell(operatorData, x, y, z);
      //       inverseCenterWeight = 1.0 / operatorData[ { 0, 0, 0 } ];
      //    }

      //    const uint_t centerIdx = vertexdof::macrocell::indexFromVertex(level, x, y, z, sd::VERTEX_C);

      //    tmp = rhs[ centerIdx ];

      //    for (const auto& neighbor : vertexdof::macrocell::neighborsWithoutCenter)
      //    {
      //       const uint_t idx = vertexdof::macrocell::indexFromVertex(level, x, y, z, neighbor);
      //       tmp -= operatorData[ vertexdof::logicalIndexOffsetFromVertex(neighbor) ] * dst[ idx ];
      //    }

      //    dst[ centerIdx ] = (1.0 - relax) * dst[ centerIdx ] + tmp * relax * inverseCenterWeight;
      // }
   }

   /// callback functions for generated kernels //////////////////////////////////////////////////

   /* apply the operator to all DoF on a given macro-face using generated kernels
   May be overridden in child-class if generated kernels are available
   */
   virtual void apply_face3D_generated( Face&                                                    face,
                                        const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcId,
                                        const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                                        const uint_t&                                            level,
                                        UpdateType                                               update ) const
   {
      apply_face3D( face, srcId, dstId, level, update );
   }

   /* apply the operator to all DoF on a given macro-face using generated kernels
      May be overridden in child-class if generated kernels are available
   */
   virtual void apply_face_generated( Face&                                                    face,
                                      const PrimitiveDataID< FunctionMemory< real_t >, Face >& srcId,
                                      const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                                      const uint_t&                                            level,
                                      UpdateType                                               update ) const
   {
      apply_face( face, srcId, dstId, level, update );
   }

   /* apply the operator to all DoF on a given macro-cell using generated kernels
      May be overridden in child-class if generated kernels are available
   */
   virtual void apply_cell_generated( Cell&                                                    cell,
                                      const PrimitiveDataID< FunctionMemory< real_t >, Cell >& srcId,
                                      const PrimitiveDataID< FunctionMemory< real_t >, Cell >& dstId,
                                      const uint_t&                                            level,
                                      UpdateType                                               update ) const
   {
      apply_cell( cell, srcId, dstId, level, update );
   }

   /* apply the SOR-operator to all DoF on a given macro-face using generated kernels
      May be overridden in child-class if generated kernels are available
   */
   virtual void smooth_sor_face3D_generated( Face&                                                    face,
                                             const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                                             const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsId,
                                             const uint_t&                                            level,
                                             real_t                                                   relax,
                                             const bool&                                              backwards = false ) const
   {
      smooth_sor_face3D( face, dstId, rhsId, level, relax, backwards );
   }

   /* apply the SOR-operator to all DoF on a given macro-face using generated kernels
      May be overridden in child-class if generated kernels are available
   */
   virtual void smooth_sor_face_generated( Face&                                                    face,
                                           const PrimitiveDataID< FunctionMemory< real_t >, Face >& dstId,
                                           const PrimitiveDataID< FunctionMemory< real_t >, Face >& rhsId,
                                           const uint_t&                                            level,
                                           real_t                                                   relax,
                                           const bool&                                              backwards = false ) const
   {
      smooth_sor_face( face, dstId, rhsId, level, relax, backwards );
   }

   /* apply the SOR-operator to all DoF on a given macro-cell using generated kernels
      May be overridden in child-class if generated kernels are available
   */
   virtual void smooth_sor_cell_generated( Cell&                                                    cell,
                                           const PrimitiveDataID< FunctionMemory< real_t >, Cell >& dstId,
                                           const PrimitiveDataID< FunctionMemory< real_t >, Cell >& rhsId,
                                           const uint_t&                                            level,
                                           real_t                                                   relax,
                                           const bool&                                              backwards = false ) const
   {
      smooth_sor_cell( cell, dstId, rhsId, level, relax, backwards );
   }

 protected:
   /// functions for variable stencil assembly. To be used in ,e.g., Ctor of constant operator, callback functions of variable operator, etc.

   /* Initialize assembly of variable edge stencil.
   */
   inline void assemble_variableStencil_edge_init( Edge& edge, const uint_t level ) const
   {
      h_           = 1.0 / ( walberla::real_c( levelinfo::num_microvertices_per_edge( level ) - 1 ) );
      level_       = level;
      stencilSize_ = edge.getData( edgeStencilID_ )->getSize( level );
      edge_        = &edge;

      // 3D version
      if ( storage_->hasGlobalCells() )
      {
         // new map not yet used -> prepare only linear stencil memory
         form_.setGeometryMap( edge.getGeometryMap() );
      }
      // 2D version
      else
      {
         x0_ = edge.getCoordinates()[0];
         dx_ = h_ * edge.getDirection();

         Face* faceS = storage_->getFace( edge.neighborFaces()[0] );
         formS_.setGeometryMap( faceS->getGeometryMap() );

         Face* faceN = nullptr;

         if ( edge.getNumNeighborFaces() == 2 )
         {
            north_ = true;
            faceN  = storage_->getFace( edge.neighborFaces()[1] );
            formN_.setGeometryMap( faceN->getGeometryMap() );
         }
         else
         {
            north_ = false;
         }

         stencil_directions_2D_ = stencil::Directions2D( h_, edge, faceS, faceN );
      }
   }

   /* assembly of variable edge stencil (requires assemble_variableStencil_edge_init() for appropriate edge and level).
   */
   inline void assemble_variableStencil_edge( real_t* edge_stencil, const uint_t i ) const
   {
      using namespace vertexdof::macroedge;
      using sD = stencilDirection;

      Point3D x = x0_ + walberla::real_t(i) * dx_;

      // 3D version (old version)
      if ( storage_->hasGlobalCells() )
      {
         // old linear stencil still used at certain points (e.g. in P2ConstantOperator)
         auto stencil = P1Elements::P1Elements3D::assembleP1LocalStencil_new< P1Form >(
             storage_, *edge_, indexing::Index( i, 0, 0 ), level_, form_ );

         WALBERLA_ASSERT_EQUAL( stencilSize_, stencil.size() );

         for ( uint_t j = 0; j < stencilSize_; j++ )
         {
            edge_stencil[j] = stencil[j];
         }

         if constexpr ( Lumped )
         {
            edge_stencil[stencilIndexOnEdge( sD::VERTEX_C )] += edge_stencil[stencilIndexOnEdge( sD::VERTEX_W )];
            edge_stencil[stencilIndexOnEdge( sD::VERTEX_W )] = 0;
            edge_stencil[stencilIndexOnEdge( sD::VERTEX_C )] += edge_stencil[stencilIndexOnEdge( sD::VERTEX_E )];
            edge_stencil[stencilIndexOnEdge( sD::VERTEX_E )] = 0;
            for ( uint_t neighborFace = 0; neighborFace < edge_->getNumNeighborFaces(); neighborFace++ )
            {
               edge_stencil[stencilIndexOnEdge( sD::VERTEX_C )] +=
                   edge_stencil[stencilIndexOnNeighborFace( sD::VERTEX_W, neighborFace )];
               edge_stencil[stencilIndexOnNeighborFace( sD::VERTEX_W, neighborFace )] = 0;
               edge_stencil[stencilIndexOnEdge( sD::VERTEX_C )] +=
                   edge_stencil[stencilIndexOnNeighborFace( sD::VERTEX_E, neighborFace )];
               edge_stencil[stencilIndexOnNeighborFace( sD::VERTEX_E, neighborFace )] = 0;
            }
            for ( uint_t neighborCell = 0; neighborCell < edge_->getNumNeighborCells(); neighborCell++ )
            {
               edge_stencil[stencilIndexOnEdge( sD::VERTEX_C )] +=
                   edge_stencil[stencilIndexOnNeighborCell( neighborCell, edge_->getNumNeighborFaces() )];
               edge_stencil[stencilIndexOnNeighborCell( neighborCell, edge_->getNumNeighborFaces() )] = 0;
            }
         }
         if constexpr ( Diagonal )
         {
            edge_stencil[stencilIndexOnEdge( sD::VERTEX_W )] = 0;
            edge_stencil[stencilIndexOnEdge( sD::VERTEX_E )] = 0;
            for ( uint_t neighborFace = 0; neighborFace < edge_->getNumNeighborFaces(); neighborFace++ )
            {
               edge_stencil[stencilIndexOnNeighborFace( sD::VERTEX_W, neighborFace )] = 0;
               edge_stencil[stencilIndexOnNeighborFace( sD::VERTEX_E, neighborFace )] = 0;
            }
            for ( uint_t neighborCell = 0; neighborCell < edge_->getNumNeighborCells(); neighborCell++ )
            {
               edge_stencil[stencilIndexOnNeighborCell( neighborCell, edge_->getNumNeighborFaces() )] = 0;
            }
         }
         if constexpr ( InvertDiagonal )
         {
            edge_stencil[stencilIndexOnEdge( sD::VERTEX_C )] = 1.0 / edge_stencil[stencilIndexOnEdge( sD::VERTEX_C )];
         }
      }
      // 2D version
      else
      {
         std::memset( edge_stencil, 0, stencilSize_ * sizeof( real_t ) );

         // south face
         vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
             formS_,
             { x, x + stencil_directions_2D_.W, x + stencil_directions_2D_.S },
             P1Elements::P1Elements2D::elementSW,
             edge_stencil );
         vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
             formS_,
             { x, x + stencil_directions_2D_.S, x + stencil_directions_2D_.SE },
             P1Elements::P1Elements2D::elementS,
             edge_stencil );
         vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
             formS_,
             { x, x + stencil_directions_2D_.SE, x + stencil_directions_2D_.E },
             P1Elements::P1Elements2D::elementSE,
             edge_stencil );

         // north face
         if ( north_ )
         {
            vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
                formN_,
                { x, x + stencil_directions_2D_.E, x + stencil_directions_2D_.N },
                P1Elements::P1Elements2D::elementNE,
                edge_stencil );
            vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
                formN_,
                { x, x + stencil_directions_2D_.N, x + stencil_directions_2D_.NW },
                P1Elements::P1Elements2D::elementN,
                edge_stencil );
            vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
                formN_,
                { x, x + stencil_directions_2D_.NW, x + stencil_directions_2D_.W },
                P1Elements::P1Elements2D::elementNW,
                edge_stencil );
         }

         if constexpr ( Lumped )
         {
            for ( const auto& neighbor : neighborsOnEdgeFromVertexDoF )
            {
               edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] +=
                   edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )];
               edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
            }

            for ( const auto& neighbor : neighborsOnSouthFaceFromVertexDoF )
            {
               edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] +=
                   edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )];
               edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
            }

            if ( north_ )
            {
               for ( const auto& neighbor : neighborsOnNorthFaceFromVertexDoF )
               {
                  edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] +=
                      edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )];
                  edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
               }
            }
         }

         if constexpr ( Diagonal )
         {
            for ( const auto& neighbor : neighborsOnEdgeFromVertexDoF )
            {
               edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
            }

            for ( const auto& neighbor : neighborsOnSouthFaceFromVertexDoF )
            {
               edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
            }

            if ( north_ )
            {
               for ( const auto& neighbor : neighborsOnNorthFaceFromVertexDoF )
               {
                  edge_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
               }
            }
         }

         if constexpr ( InvertDiagonal )
         {
            edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] =
                1.0 / edge_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )];
         }
      }
   }

   /* assembly of variable edge stencil (requires assemble_variableStencil_edge_init() for appropriate edge and level).
   */
   inline void assemble_variableStencil_edge3D( vertexdof::macroedge::StencilMap_T& edge_stencil, const uint_t i ) const
   {
      WALBERLA_ASSERT( storage_->hasGlobalCells() );

      for ( uint_t neighborCellID = 0; neighborCellID < edge_->getNumNeighborCells(); ++neighborCellID )
      {
         auto neighborCell = storage_->getCell( edge_->neighborCells().at( neighborCellID ) );
         auto vertexAssemblyIndexInCell =
             vertexdof::macroedge::getIndexInNeighboringMacroCell( { i, 0, 0 }, *edge_, neighborCellID, *storage_, level_ );
         edge_stencil[neighborCellID] = P1Elements::P1Elements3D::assembleP1LocalStencilNew_new< P1Form >(
             storage_, *neighborCell, vertexAssemblyIndexInCell, level_, form_ );
      }

      // todo: also implement lumping/diagonal/invertdiagonal here
      /* remark: This has also been missing in the old implementation of P1ConstantOperator.
            Since assemble_variableStencil_edge3D() is only used in P2ConstantOperator,
            it is not required for now.
      */
   }

   /* Initialize assembly of variable face stencil.
   */
   inline void assemble_variableStencil_face_init( Face& face, const uint_t level ) const
   {
      h_           = 1.0 / ( walberla::real_c( levelinfo::num_microvertices_per_edge( level ) - 1 ) );
      level_       = level;
      stencilSize_ = face.getData( faceStencilID_ )->getSize( level );
      face_        = &face;

      form_.setGeometryMap( face.getGeometryMap() );

      // 3D version
      if ( storage_->hasGlobalCells() )
      {
         // nothing to do here
      }
      // 2D version
      else
      {
         x0_ = face.coords[0];
         dx_ = h_ * ( face.coords[1] - face.coords[0] );
         dy_ = h_ * ( face.coords[2] - face.coords[0] );

         stencil_directions_2D_ = stencil::Directions2D( h_, face );
      }
   }

   /* assembly of variable face stencil (requires assemble_variableStencil_face_init() for appropriate face and level).
   */
   inline void assemble_variableStencil_face( real_t* face_stencil, const uint_t i, const uint_t j ) const
   {
      using sD = stencilDirection;

      std::memset( face_stencil, 0, stencilSize_ * sizeof( real_t ) );

      WALBERLA_ASSERT( !storage_->hasGlobalCells() );
      Point3D x = x0_ + walberla::real_t(i) * dx_ + walberla::real_t(j) * dy_;

      vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
          form_,
          { x, x + stencil_directions_2D_.W, x + stencil_directions_2D_.S },
          P1Elements::P1Elements2D::elementSW,
          face_stencil );
      vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
          form_,
          { x, x + stencil_directions_2D_.S, x + stencil_directions_2D_.SE },
          P1Elements::P1Elements2D::elementS,
          face_stencil );
      vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
          form_,
          { x, x + stencil_directions_2D_.SE, x + stencil_directions_2D_.E },
          P1Elements::P1Elements2D::elementSE,
          face_stencil );
      vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
          form_,
          { x, x + stencil_directions_2D_.E, x + stencil_directions_2D_.N },
          P1Elements::P1Elements2D::elementNE,
          face_stencil );
      vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
          form_,
          { x, x + stencil_directions_2D_.N, x + stencil_directions_2D_.NW },
          P1Elements::P1Elements2D::elementN,
          face_stencil );
      vertexdof::variablestencil::assembleLocalStencil_new< P1Form >(
          form_,
          { x, x + stencil_directions_2D_.NW, x + stencil_directions_2D_.W },
          P1Elements::P1Elements2D::elementNW,
          face_stencil );

      if constexpr ( Lumped )
      {
         for ( const auto& neighbor : vertexdof::macroface::neighborsWithoutCenter )
         {
            face_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] +=
                face_stencil[vertexdof::stencilIndexFromVertex( neighbor )];
            face_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
         }
      }

      if constexpr ( Diagonal )
      {
         for ( const auto& neighbor : vertexdof::macroface::neighborsWithoutCenter )
         {
            face_stencil[vertexdof::stencilIndexFromVertex( neighbor )] = 0;
         }
      }

      if constexpr ( InvertDiagonal )
      {
         face_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )] =
             1.0 / face_stencil[vertexdof::stencilIndexFromVertex( sD::VERTEX_C )];
      }
   }

   /* assembly of variable face stencil (requires assemble_variableStencil_face_init() for appropriate face and level).
   */
   inline void
       assemble_variableStencil_face3D( vertexdof::macroface::StencilMap_T& face_stencil, const uint_t i, const uint_t j ) const
   {
      WALBERLA_ASSERT( storage_->hasGlobalCells() );

      for ( uint_t neighborCellID = 0; neighborCellID < face_->getNumNeighborCells(); ++neighborCellID )
      {
         auto neighborCell = storage_->getCell( face_->neighborCells().at( neighborCellID ) );
         auto vertexAssemblyIndexInCell =
             vertexdof::macroface::getIndexInNeighboringMacroCell( { i, j, 0 }, *face_, neighborCellID, *storage_, level_ );
         face_stencil[neighborCellID] = P1Elements::P1Elements3D::assembleP1LocalStencilNew_new< P1Form >(
             storage_, *neighborCell, vertexAssemblyIndexInCell, level_, form_ );
      }

      if constexpr ( Lumped )
      {
         for ( uint_t neighborCellID = 0; neighborCellID < face_->getNumNeighborCells(); neighborCellID++ )
         {
            for ( auto& stencilIt : face_stencil[neighborCellID] )
            {
               if ( !( neighborCellID == 0 && stencilIt.first == indexing::IndexIncrement( { 0, 0, 0 } ) ) )
               {
                  face_stencil[0][{ 0, 0, 0 }] += stencilIt.second;
                  stencilIt.second = 0;
               }
            }
         }
      }
      if constexpr ( Diagonal )
      {
         for ( uint_t neighborCellID = 0; neighborCellID < face_->getNumNeighborCells(); neighborCellID++ )
         {
            for ( auto& stencilIt : face_stencil[neighborCellID] )
            {
               if ( stencilIt.first != indexing::IndexIncrement( { 0, 0, 0 } ) )
               {
                  stencilIt.second = 0;
               }
            }
         }
      }
      if constexpr ( InvertDiagonal )
      {
         for ( uint_t neighborCellID = 1; neighborCellID < face_->getNumNeighborCells(); neighborCellID++ )
         {
            face_stencil[0][{ 0, 0, 0 }] += face_stencil[neighborCellID][{ 0, 0, 0 }];
            face_stencil[neighborCellID][{ 0, 0, 0 }] = 0;
         }
         face_stencil[0][{ 0, 0, 0 }] = 1.0 / face_stencil[0][{ 0, 0, 0 }];
      }
   }

   /* Initialize assembly of variable cell stencil.
   */
   inline void assemble_variableStencil_cell_init( Cell& cell, const uint_t level ) const
   {
      h_ = 1.0 / ( walberla::real_c( levelinfo::num_microvertices_per_edge( level ) - 1 ) );
      form_.setGeometryMap( cell.getGeometryMap() );
      level_ = level;
      cell_  = &cell;
   }

   /* assembly of variable cell stencil (requires assemble_variableStencil_cell_init() for appropriate cell and level).
   */
   inline void assemble_variableStencil_cell( vertexdof::macrocell::StencilMap_T& cell_stencil,
                                              const uint_t                        i,
                                              const uint_t                        j,
                                              const uint_t                        k ) const
   {
      cell_stencil = P1Elements::P1Elements3D::assembleP1LocalStencilNew_new< P1Form >(
          storage_, *cell_, indexing::Index( i, j, k ), level_, form_ );

      if constexpr ( Lumped )
      {
         for ( auto dir : vertexdof::macrocell::neighborsWithoutCenter )
         {
            cell_stencil[{ 0, 0, 0 }] += cell_stencil[vertexdof::logicalIndexOffsetFromVertex( dir )];
            cell_stencil[vertexdof::logicalIndexOffsetFromVertex( dir )] = 0;
         }
      }
      if constexpr ( Diagonal )
      {
         for ( auto dir : vertexdof::macrocell::neighborsWithoutCenter )
         {
            cell_stencil[vertexdof::logicalIndexOffsetFromVertex( dir )] = 0;
         }
      }
      if constexpr ( InvertDiagonal )
      {
         cell_stencil[{ 0, 0, 0 }] = 1.0 / cell_stencil[{ 0, 0, 0 }];
      }
   }

   ////////////////////////////////////////////////////////////////////////////////////////////////////////

   // return true if backwards sor is implemented for this operator, false otherwise
   virtual bool backwards_sor_available() const = 0;

   // return true if the stencil has to be recomputed for every DoF
   virtual bool variableStencil() const = 0;

   std::shared_ptr< P1Function< real_t > > diagonalValues_;
   std::shared_ptr< P1Function< real_t > > inverseDiagonalValues_;

   PrimitiveDataID< StencilMemory< real_t >, Vertex >                             vertexStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Edge >                               edgeStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Face >                               faceStencilID_;
   PrimitiveDataID< LevelWiseMemory< vertexdof::macroedge::StencilMap_T >, Edge > edgeStencil3DID_;
   PrimitiveDataID< LevelWiseMemory< vertexdof::macroface::StencilMap_T >, Face > faceStencil3DID_;
   PrimitiveDataID< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell > cellStencilID_;

   // general data for stencil assembly
   mutable Point3D x0_, dx_, dy_, dz_;
   mutable uint_t  level_;
   mutable real_t  h_;
   mutable P1Form  form_;

   // data for edge stencil assembly
   mutable stencil::Directions2D stencil_directions_2D_;
   mutable P1Form                formS_, formN_;
   mutable bool                  north_;
   mutable Edge*                 edge_;
   mutable uint_t                stencilSize_;

   // data for face stencil assembly
   mutable Face* face_;

   // data for cell stencil assembly
   mutable Cell* cell_;

   /// VIRTUAL CALLBACK FUNCTIONS -- TO BE IMPLEMENTED IN CHILD CLASSES (e.g. constant, variable, surrogate, ...) ///

   /* Initialize assembly of variable edge stencil.
      Will be called before iterating over edge whenever the stencil is applied.
   */
   virtual void assemble_stencil_edge_init( Edge& edge, const uint_t level ) const = 0;

   /* Assembly of edge stencil.
      Will be called before stencil is applied to a particuar edge-DoF.
   */
   virtual void assemble_stencil_edge( real_t* edge_stencil, const uint_t i ) const = 0;

   /* Initialize assembly of face stencil.
      Will be called before iterating over face whenever the stencil is applied.
   */
   virtual void assemble_stencil_face_init( Face& face, const uint_t level ) const = 0;

   /* Initialize data before continuing iteration over DoF in x direction
      (required for surrogate polynomials)
   */
   virtual void assemble_stencil_face_init_y( const uint_t j ) const {};

   /* Assembly of face stencil.
      Will be called before stencil is applied to a particuar face-DoF of a 2d domain.
   */
   virtual void assemble_stencil_face( real_t* face_stencil, const uint_t i, const uint_t j ) const = 0;

   /* Assembly of face stencil.
      Will be called before stencil is applied to a particuar face-DoF of a 3D domain.
   */
   virtual void
       assemble_stencil_face3D( vertexdof::macroface::StencilMap_T& face_stencil, const uint_t i, const uint_t j ) const = 0;

   /* Initialize assembly of cell stencil.
      Will be called before iterating over cell whenever the stencil is applied.
   */
   virtual void assemble_stencil_cell_init( Cell& cell, const uint_t level ) const = 0;

   /* Initialize data before continuing iteration over DoF in y direction
      (required for surrogate polynomials)
   */
   virtual void assemble_stencil_cell_init_z( const uint_t k ) const {};

   /* Initialize data before continuing iteration over DoF in x direction
      (required for surrogate polynomials)
   */
   virtual void assemble_stencil_cell_init_y( const uint_t j ) const {};

   /* Assembly of cell stencil.
      Will be called before stencil is applied to a particuar cell-DoF.
   */
   virtual void assemble_stencil_cell( vertexdof::macrocell::StencilMap_T& cell_stencil,
                                       const uint_t                        i,
                                       const uint_t                        j,
                                       const uint_t                        k ) const = 0;
};

} // namespace hyteg
