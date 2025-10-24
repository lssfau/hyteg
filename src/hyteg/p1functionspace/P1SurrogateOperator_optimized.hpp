/*
 * Copyright (c) 2025 Benjamin Mann
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
#include <hyteg/Stencil.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q2.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_affine_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_blending_q3.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_k_mass_affine_q4.hpp>
#include <hyteg/forms/form_hyteg_generated/p1/p1_mass_blending_q4.hpp>
#include <hyteg/operators/Operator.hpp>
#include <hyteg/p1functionspace/P1Elements.hpp>
#include <hyteg/p1functionspace/P1Function.hpp>
#include <hyteg/p1functionspace/VertexDoFMacroFace.hpp>
#include <hyteg/p1functionspace/globalIndices.hpp>
#include <hyteg/polynomial/elementwise/data.hpp>
#include <hyteg/polynomial/stencil/leastSquares.hpp>
#include <hyteg/polynomial/stencil/polynomial.hpp>
#include <hyteg/solvers/Smoothables.hpp>

namespace hyteg {

template < class P1Form, uint8_t DEGREE >
class P1SurrogateOperator : public Operator< P1Function< real_t >, P1Function< real_t > >,
                            public GSSmoothable< P1Function< real_t > >,
                            public GSBackwardsSmoothable< P1Function< real_t > >,
                            public SORSmoothable< P1Function< real_t > >,
                            public SORBackwardsSmoothable< P1Function< real_t > >,
                            public WeightedJacobiSmoothable< P1Function< real_t > >,
                            public OperatorWithInverseDiagonal< P1Function< real_t > >
{
   /* On lower levels, storing and evaluating polynomials is significantly less performant.
      Therefore, on levels 0-3 we precompute and store the system matrices, while we use
      surrogates for levels 4+
    */
   static constexpr uint_t min_lvl_for_surrogate = 4;

   /* Single precision LSQ may lead to very poor accuracy of the resulting polynomials. Therefore,
      we use real_t for the polynomial evaluation only, while sticking to double precision LSQ.
    */
   using PolyDomain = surrogate::polynomial::Domain< real_t >;
   using LSQ        = p1::stencil::surrogate::LeastSquares< double >;

   // regular stencils, i.e., cell-dof (3d), face-dof (2d/3d), edge-dof (2d)
   template < uint8_t DIM >
   using Stencil = p1::stencil::StencilData< DIM >;
   // surrogate stencils
   template < uint8_t DIM_domain, uint8_t DIM_primitive >
   using PolyStencil = p1::stencil::surrogate::Polynomial< DIM_domain, DIM_primitive, DEGREE >;
   // irregular stencils, i.e., edge-dof (3d), vertex-dof (2d/3d)
   using VarStencil = std::vector< real_t >;
   // stencil data for LSQ fit
   template < uint8_t DIM >
   using LSQData = p1::stencil::StencilData< DIM, typename surrogate::LeastSquares< real_t >::Vector >;
   // containers for precomputed stencils. usage: map[id][lvl][k] -> stencil at k-th DoF of of el_id on lvl
   template < uint8_t DIM >
   using StencilMap    = surrogate::ElementWiseData< std::vector< Stencil< DIM > > >;
   using VarStencilMap = surrogate::ElementWiseData< std::vector< VarStencil > >;
   // container for surrogate stencils. usage: map[id][lvl] -> poly-stencil approximating the stencil of el_id on lvl
   template < uint8_t DIM_domain, uint8_t DIM_primitive >
   using PolyStencilMap = surrogate::ElementWiseData< PolyStencil< DIM_domain, DIM_primitive > >;
   // data structure for DoF indices
   using DofIdx = p1::stencil::StencilData< 3, walberla::uint_t >;

 public:
   P1SurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : P1SurrogateOperator( storage, minLevel, maxLevel, P1Form() )
   {}

   P1SurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel, const P1Form& form )
   : Operator< P1Function< real_t >, P1Function< real_t > >( storage, minLevel, maxLevel )
   , form_( form )
   , is_initialized_( false )
   , lsq_volume_( maxLevel + 1 )
   , lsq_interface_( maxLevel + 1 )
   , downsampling_( maxLevel + 1 )
   , stencil_vtx_( storage, maxLevel, 0 )
   , stencil_edge_3d_( storage, maxLevel, 1 )
   , stencil_edge_2d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ), 1 )
   , stencil_face_2d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ), 2 )
   , stencil_face_3d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ), 2 )
   , stencil_cell_3d_( storage, std::min( maxLevel, min_lvl_for_surrogate - 1u ), 3 )
   , surrogate_edge_2d_( storage, maxLevel, 1 )
   , surrogate_face_2d_( storage, maxLevel, 2 )
   , surrogate_face_3d_( storage, maxLevel, 2 )
   , surrogate_cell_3d_( storage, maxLevel, 3 )
   {}

   init( size_t downsampling, const std::string& path_to_svd, bool needsInverseDiagEntries )
   {
      uint_t dim = ( storage_->hasGlobalCells() ) ? 3 : 2;

      if ( dim == 3 )
      {
         // mapping between logical indices and stencil directions for each face
         compute_offsets_for_face_stencil_3d();
      }

      /* precompute and store stencils
         * irregular stencils are precomputed for all levels
         * regular stencils are precomputed for levels 1-3
      */
      for ( uint_t level = 0; level <= maxLevel_; ++level )
      {
         if ( dim == 2 )
         {
            precompute_stencil_vtx_2d( level );
            if ( 1 <= level && level < min_lvl_for_surrogate )
            {
               precompute_stencil_edge_2d( level );
            }
            if ( 2 <= level && level < min_lvl_for_surrogate )
            {
               precompute_stencil_face_2d( level );
            }
         }
         else
         {
            precompute_stencil_vtx_3d( level );
            if ( 1 <= level )
            {
               precompute_stencil_edge_3d( level );
            }
            if ( 2 <= level && level < min_lvl_for_surrogate )
            {
               precompute_stencil_face_3d( level );
               precompute_stencil_cell_3d( level );
            }
         }
      }

      // approximate regular stencils for level 4+ by polynomials
      for ( uint_t level = min_lvl_for_surrogate; level <= maxLevel_; ++level )
      {
         // initialize least squares approximation
         if ( lsq_volume_[level] == nullptr || downsampling_ != downsampling )
         {
            if ( path_to_svd == "" )
            {
               lsq_volume_[level]    = std::make_shared< LSQ >( dim, DEGREE, level, downsampling );
               lsq_interface_[level] = std::make_shared< LSQ >( dim - 1, DEGREE, level, downsampling );
            }
            else
            {
               lsq_volume_[level]    = std::make_shared< LSQ >( path_to_svd, dim, DEGREE, level, downsampling );
               lsq_interface_[level] = std::make_shared< LSQ >( path_to_svd, dim - 1, DEGREE, level, downsampling );
            }
            downsampling_ = downsampling;
         }

         if ( dim == 2 )
         {
            compute_surrogates_edge_2d( level );
            compute_surrogates_face_2d( level );
         }
         else
         {
            compute_surrogates_face_3d( level );
            compute_surrogates_cell_3d( level );
         }
      }

      if ( needsInverseDiagEntries )
      {
         computeInverseDiagonalOperatorValues();
      }

      is_initialized_ = true;
   }

   /* compute h^(d/2)*||A - Aq||_F with variable operator A and surrogate Aq
      @returns [h^(d/2)*||A - Aq||_F restricted to K] for all macro elements K
   */
   std::map< PrimitiveID, real_t > computeSurrogateError( uint_t level ) const
   {
      if ( storage_->hasGlobalCells() )
      {
         return computeSurrogateError3D( level );
      }
      else
      {
         return computeSurrogateError2D( level );
      }
   }

   store_svd( const std::string& path_to_svd )
   {
      for ( uint_t level = min_lvl_for_surrogate; level <= maxLevel_; ++level )
      {
         if ( lsq_volume_[level] != nullptr )
         {
            lsq_volume_[level]->write_to_file( path_to_svd );
            lsq_interface_[level]->write_to_file( path_to_svd );
         }
      }
   }

   void apply( const P1Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const override final
   {
      WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

      this->startTiming( "Apply" );
      src.template communicate< Vertex, Edge >( level );
      src.template communicate< Edge, Face >( level );
      src.template communicate< Face, Cell >( level );

      src.template communicate< Cell, Face >( level );
      src.template communicate< Face, Edge >( level );
      src.template communicate< Edge, Vertex >( level );

      const int dim = ( storage_->hasGlobalCells() ) ? 3 : 2;

      this->timingTree_->start( "Macro-Vertex" );

      for ( const auto& [vtxId, vtx] : storage_->getVertices() )
      {
         const auto bc = dst.getBoundaryCondition().getBoundaryType( vtx->getMeshBoundaryFlag() );

         if ( testFlag( bc, flag ) )
         {
            auto        srcData = vtx->getData( src.getVertexDataID() )->getPointer( level );
            auto        dstData = vtx->getData( dst.getVertexDataID() )->getPointer( level );
            const auto& stencil = stencil_vtx_.at( vtxId )[level][0];

            if ( updateType == Replace )
            {
               dstData[0] = real_c( 0 );
            } // else updateType == Add

            for ( uint_t i = 0; i < stencil.size(); ++i )
            {
               dstData[0] += stencil[i] * srcData[i];
            }
         }
      }

      this->timingTree_->stop( "Macro-Vertex" );

      this->timingTree_->start( "Macro-Edge" );

      if ( level >= 1 )
      {
         for ( const auto& [edgeId, edge] : storage_->getEdges() )
         {
            const auto bc = dst.getBoundaryCondition().getBoundaryType( edge->getMeshBoundaryFlag() );

            if ( testFlag( bc, flag ) )
            {
               auto srcData = edge->getData( src.getEdgeDataID() )->getPointer( level );
               auto dstData = edge->getData( dst.getEdgeDataID() )->getPointer( level );

               if ( dim == 2 )
               {
                  if ( level < min_lvl_for_surrogate )
                  {
                     apply_edge_precomputed_2d( edge, level, srcData, dstData );
                  }
                  else
                  {
                     apply_edge_surrogate_2d( edge, level, srcData, dstData );
                  }
               }
               else // dim == 3
               {
                  apply_edge_precomputed_3d( edge, level, srcData, dstData );
               }
            }
         }
      }

      this->timingTree_->stop( "Macro-Edge" );

      this->timingTree_->start( "Macro-Face" );

      if ( level >= 2 )
      {
         for ( const auto& [faceId, face] : storage_->getFaces() )
         {
            const auto bc = dst.getBoundaryCondition().getBoundaryType( face->getMeshBoundaryFlag() );

            if ( testFlag( bc, flag ) )
            {
               auto srcData = face->getData( src.getFaceDataID() )->getPointer( level );
               auto dstData = face->getData( dst.getFaceDataID() )->getPointer( level );

               if ( dim == 2 )
               {
                  if ( level < min_lvl_for_surrogate )
                  {
                     apply_face_precomputed_2d( face, level, srcData, dstData );
                  }
                  else
                  {
                     apply_face_surrogate_2d( face, level, srcData, dstData );
                  }
               }
               else // dim == 3
               {
                  if ( level < min_lvl_for_surrogate )
                  {
                     apply_face_precomputed_3d( face, level, srcData, dstData );
                  }
                  else
                  {
                     apply_face_surrogate_3d( face, level, srcData, dstData );
                  }
               }
            }
         }
      }

      this->timingTree_->stop( "Macro-Face" );

      this->timingTree_->start( "Macro-Cell" );

      if ( level >= 2 )
      {
         for ( const auto& [cellId, cell] : storage_->getCells() )
         {
            const auto bc = dst.getBoundaryCondition().getBoundaryType( cell->getMeshBoundaryFlag() );

            if ( testFlag( bc, flag ) )
            {
               auto srcData = cell->getData( src.getCellDataID() )->getPointer( level );
               auto dstData = cell->getData( dst.getCellDataID() )->getPointer( level );

               if ( level < min_lvl_for_surrogate )
               {
                  apply_cell_precomputed_3d( cell, level, srcData, dstData );
               }
               else
               {
                  apply_cell_surrogate_3d( cell, level, srcData, dstData );
               }
            }
         }
      }

      this->timingTree_->stop( "Macro-Cell" );

      this->stopTiming( "Apply" );
   }

   void smooth_gs( const P1Function< real_t >& dst,
                   const P1Function< real_t >& rhs,
                   size_t                      level,
                   DoFType                     flag ) const override final

   {
      smooth_sor( dst, rhs, 1.0, level, flag );
   }

   void smooth_sor( const P1Function< real_t >& dst,
                    const P1Function< real_t >& rhs,
                    real_t                      relax,
                    size_t                      level,
                    DoFType                     flag ) const override
   {
      // todo
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
      dst.assign( { real_t( 1.0 ), real_t( -1.0 ) }, { rhs, dst }, level, flag );

      // perform Jacobi update step
      dst.multElementwise( { *getInverseDiagonalValues(), dst }, level, flag );
      dst.assign( { static_cast< real_t >( 1.0 ), relax }, { src, dst }, level, flag );

      this->stopTiming( "smooth_jac" );
   }

   void computeDiagonalOperatorValues() { computeDiagonalOperatorValues( false ); }

   void computeInverseDiagonalOperatorValues() override final { computeDiagonalOperatorValues( true ); }

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

 private:
   std::map< PrimitiveID, real_t > computeSurrogateError2D( uint_t level ) const
   {
      // todo
   }

   std::map< PrimitiveID, real_t > computeSurrogateError3D( uint_t level ) const
   {
      // todo
   }

   void precompute_stencil_vtx_2d( uint_t lvl )
   {
      const auto n = levelinfo::num_microvertices_per_edge( lvl );
      const auto h = real_t( 1.0 / ( real_t( n - 1 ) ) );

      for ( const auto& [vtxId, vtx] : storage_->getVertices() )
      {
         auto& stencil = stencil_vtx_[vtxId][lvl][0];
         stencil.resize( vtx->getNumNeighborEdges() + 1 );
         std::fill( stencil.begin(), stencil.end(), real_t( 0 ) );

         for ( auto& faceId : vtx->neighborFaces() )
         {
            Face* face = storage_->getFace( faceId );
            form_.setGeometryMap( face->getGeometryMap() );

            // local vertex indices
            const auto  vx        = face->vertex_index( vtxId );
            const auto  adj_edges = face->adjacent_edges( vtxId );
            const auto  v0        = face->vertex_index( storage_->getEdge( adj_edges[0] )->get_opposite_vertex( vtxId ) )]
            const auto  v1        = face->vertex_index( storage_->getEdge( adj_edges[1] )->get_opposite_vertex( vtxId ) )]

            // compute local stiffness matrix
            const auto& faceCoords = face->getCoordinates();
            Matrixr< 1, 3 > matrixRow;
            const auto&     x  = faceCoords[vx];
            const auto      d0 = h * ( faceCoords[v0] - x );
            const auto      d1 = h * ( faceCoords[v1] - x );
            form_.integrateRow( 0, { { x, x + d0, x + d1 } }, matrixRow );

            // stencil directions
            const uint_t ix = 0;
            const uint_t i0 = vtx->edge_index( adj_edges[0] ) + 1;
            const uint_t i1 = vtx->edge_index( adj_edges[1] ) + 1;

            // add contributions from this face to stencil
            stencil[ix] += real_t( matrixRow( 0, 0 ) );
            stencil[i0] += real_t( matrixRow( 0, 1 ) );
            stencil[i1] += real_t( matrixRow( 0, 2 ) );
         }
      }
   }

   void precompute_stencil_vtx_3d( uint_t lvl )
   {
      const auto n = levelinfo::num_microvertices_per_edge( lvl );
      const auto h = real_t( 1.0 / ( real_t( n - 1 ) ) );

      for ( const auto& [vtxId, vtx] : storage_->getVertices() )
      {
         auto& stencil = stencil_vtx_[vtxId][lvl][0];
         stencil.resize( vtx->getNumNeighborEdges() + 1 );
         std::fill( stencil.begin(), stencil.end(), real_t( 0 ) );

         // neighboring vertices and edges
         auto                       x          = vtx->getCoordinates();
         auto                       nbrEdgeIds = vtx->neighborEdges();
         std::vector< PrimitiveID > nbrVtxIds( stencil.size() - 1 );
         std::vector< Point3D >     nbrMicroCoords( stencil.size() - 1 );
         for ( uint_t d = 1; d < stencil.size(); ++d )
         {
            Edge* edge            = storage_->getEdge( nbrEdgeIds[d - 1] );
            nbrVtxIds[d - 1]      = edge->get_opposite_vertex( vtxId );
            auto nbr_x            = storage_->getVertex( nbrVtxIds[d - 1] )->getCoordinates();
            nbrMicroCoords[d - 1] = h * nbr_x + ( 1 - h ) * x;
         }

         for ( auto& cellId : vtx->neighborCells() )
         {
            Cell* cell = storage_->getCell( cellId );

            // stencil directions associated with cell vertices
            std::array< uint_t, 3 > localStencilContributions;
            uint_t                  dir = 1;
            for ( auto& cellNbrVtxId : cell->neighborVertices() )
            {
               if ( cellNbrVtxId == vtxId )
               {
                  continue;
               }

               for ( uint_t j = 0; j < nbrVtxIds.size(); ++j )
               {
                  if ( cellNbrVtxId == nbrVtxIds[j] )
                  {
                     localStencilContributions[dir] = j;
                     break;
                  }
               }

               ++dir;
               if ( dir == 3 )
               {
                  break;
               }
            }

            // compute local stiffness matrix
            Matrixr< 1, 4 >          matrixRow;
            std::array< Point3D, 4 > coords{ x,
                                             nbrMicroCoords[localStencilContributions[0]],
                                             nbrMicroCoords[localStencilContributions[1]],
                                             nbrMicroCoords[localStencilContributions[2]] };
            form_.setGeometryMap( cell->getGeometryMap() );
            form_.integrateRow( 0, coords, matrixRow );

            // add contributions from this face to stencil
            stencil[0] += real_t( matrixRow( 0, 0 ) );
            for ( uint_t dir = 0; dir < 3; ++dir )
            {
               stencil[localStencilContributions[dir]] += real_t( matrixRow( 0, dir + 1 ) );
            }
         }
      }
   }

   void precompute_stencil_edge_2d( uint_t lvl )
   {
      const auto n = levelinfo::num_microvertices_per_edge( lvl );
      const auto h = real_t( 1.0 / ( real_t( n - 1 ) ) );

      P1Form      form_N( form_ );
      Face*       face_S;
      Face*       face_N;
      PrimitiveID vtxId_N, vtxId_S;

      for ( const auto& [edgeId, edge] : storage_->getEdges() )
      {
         auto& stencils = stencil_edge_2d_[edgeId][lvl];
         stencils.resize( n - 2 );

         // neighbor face
         auto n_nbr_faces = edge->getNumNeighborFaces();
         face_S           = storage_->getFace( edge->neighborFaces()[0] );
         form_.setGeometryMap( face_S->getGeometryMap() );
         vtxId_S = face_S->get_vertex_opposite_to_edge( edgeId );
         if ( n_nbr_faces == 2 )
         {
            face_N = storage_->getFace( edge->neighborFaces()[1] );
            form_N.setGeometryMap( face_N->getGeometryMap() );
            vtxId_N = face_N->get_vertex_opposite_to_edge( edgeId );
         }

         // coordinates
         auto x0 = edge->getCoordinates()[0];
         auto x1 = edge->getCoordinates()[1];
         auto dx = h * edge->getDirection();
         // coordinate offsets of stencil nbrs
         p1::stencil::StencilData< 2, Point3D > d;
         d[p1::stencil::W] = -dx;
         d[p1::stencil::E] = dx;
         if ( n_nbr_faces == 2 )
         {
            auto coord_N      = storage_->getVertex( vtxId_N )->getCoordinates();
            d[p1::stencil::N] = h * ( coord_N - x0 );
         }
         auto coord_S       = storage_->getVertex( vtxId_S )->getCoordinates();
         d[p1::stencil::S]  = h * ( coord_S - x1 );
         d[p1::stencil::NW] = d[p1::stencil::N] + d[p1::stencil::W];
         d[p1::stencil::SE] = d[p1::stencil::S] + d[p1::stencil::E];

         // loop over inner vertices on the macro edge
         for ( uint_t i = 1; i < n - 1; ++i )
         {
            auto x = x0 + real_t( i ) * dx;

            // compute local stiffness matrices and add contributions to stencil
            Matrixr< 1, 3 > matrixRow;
            auto&           stencil = stencils[i - 1];
            form_.integrateRow( 0, { { x, x + d[p1::stencil::W], x + d[p1::stencil::S] } }, matrixRow );
            stencil[p1::stencil::C] = real_t( matrixRow( 0, 0 ) );
            stencil[p1::stencil::W] = real_t( matrixRow( 0, 1 ) );
            stencil[p1::stencil::S] = real_t( matrixRow( 0, 2 ) );
            form_.integrateRow( 0, { { x, x + d[p1::stencil::S], x + d[p1::stencil::SE] } }, matrixRow );
            stencil[p1::stencil::C] += real_t( matrixRow( 0, 0 ) );
            stencil[p1::stencil::S] += real_t( matrixRow( 0, 1 ) );
            stencil[p1::stencil::SE] = real_t( matrixRow( 0, 2 ) );
            form_.integrateRow( 0, { { x, x + d[p1::stencil::SE], x + d[p1::stencil::E] } }, matrixRow );
            stencil[p1::stencil::C] += real_t( matrixRow( 0, 0 ) );
            stencil[p1::stencil::SE] += real_t( matrixRow( 0, 1 ) );
            stencil[p1::stencil::E] = real_t( matrixRow( 0, 2 ) );
            if ( n_nbr_faces == 2 )
            {
               form_N.integrateRow( 0, { { x, x + d[p1::stencil::E], x + d[p1::stencil::N] } }, matrixRow );
               stencil[p1::stencil::C] += real_t( matrixRow( 0, 0 ) );
               stencil[p1::stencil::E] += real_t( matrixRow( 0, 1 ) );
               stencil[p1::stencil::N] = real_t( matrixRow( 0, 2 ) );
               form_N.integrateRow( 0, { { x, x + d[p1::stencil::N], x + d[p1::stencil::NW] } }, matrixRow );
               stencil[p1::stencil::C] += real_t( matrixRow( 0, 0 ) );
               stencil[p1::stencil::N] += real_t( matrixRow( 0, 1 ) );
               stencil[p1::stencil::NW] = real_t( matrixRow( 0, 2 ) );
               form_N.integrateRow( 0, { { x, x + d[p1::stencil::NW], x + d[p1::stencil::W] } }, matrixRow );
               stencil[p1::stencil::C] += real_t( matrixRow( 0, 0 ) );
               stencil[p1::stencil::NW] += real_t( matrixRow( 0, 1 ) );
               stencil[p1::stencil::W] += real_t( matrixRow( 0, 2 ) );
            }
         }
      }
   }

   void precompute_stencil_edge_3d( uint_t lvl )
   {
      const auto n = levelinfo::num_microvertices_per_edge( lvl );
      const auto h = real_t( 1.0 / ( real_t( n - 1 ) ) );

      const indexing::Index idx1( 1, 0, 0 );
      Matrixr< 1, 4 >       localStiffnessMatrixRow;

      for ( const auto& [edgeId, edge] : storage_->getEdges() )
      {
         auto& stencils = stencil_edge_3d_[edgeId][lvl];
         stencils.resize( n - 2 );
         const auto stencil_size = vertexDoFMacroEdgeStencilMemorySize( lvl, edge );
         const auto n_nbr_faces  = edge->getNumNeighborFaces();
         const auto n_nbr_cells  = edge->getNumNeighborCells();

         // collect required data for each neighbor cell
         const auto&                                            nbrCells = edge->neighborCells();
         std::vector< P1Form >                                  form( n_nbr_cells, form_ );
         std::vector< uint_t >                                  n_microCells( n_nbr_cells );
         std::vector< std::vector< std::array< Point3D, 4 > > > coordinateOffset( n_nbr_cells );
         std::vector< std::vector< std::array< uint_t, 4 > > >  stencilIndex( n_nbr_cells );

         for ( uint_t c = 0; c < n_nbr_cells; ++c )
         {
            const Cell* cell = storage_->getCell( nbrCells[c] );
            form[c].setGeometryMap( cell->getGeometryMap() );

            // local coordinate system
            const auto              localEdgeId = cell->getLocalEdgeID( edgeId );
            std::array< uint_t, 4 > indexingBasis{ 0, 0, 5, 5 };
            // origin
            indexingBasis[0] = cell->getEdgeLocalVertexToCellLocalVertexMaps()[localEdgeID].at( 0 );
            // x-direction
            indexingBasis[1] = cell->getEdgeLocalVertexToCellLocalVertexMaps()[localEdgeID].at( 1 );
            // find y and z direction
            for ( uint_t d = 0; d < 4; ++d )
            {
               if ( d != indexingBasis[0] && d != indexingBasis[1] )
               {
                  if ( indexingBasis[2] > 4 )
                     indexingBasis[2] = d;
                  else if ( indexingBasis[3] > 4 )
                     indexingBasis[3] = d;
               }
            }
            const auto idx1_cell = indexing::basisConversion( idx1, indexingBasis, { 0, 1, 2, 3 }, n );
            const auto center    = vertexdof::macrocell::coordinateFromIndex( lvl, *cell, idx1_cell );

            // intersection of nbrFaces from edge and cell
            auto edgeFaces = vertexdof::macrocell::isOnCellFace( idx1_cell, lvl );

            // micro-cells associated with the stencil
            const auto microCells = P1Elements::P1Elements3D::getNeighboringElements( idx1_cell, lvl );
            n_microCells[c]       = microCells.size();
            coordinateOffset[c].resize( microCells.size() );
            stencilIndex[c].resize( microCells.size() );

            // iterate over micro-cells
            for ( uint_t mc = 0; mc < microCells.size(); ++mc )
            {
               const auto& microCell = microCells[mc];
               // iterate over the micro-vertices of the micro-cell
               for ( uint_t mv = 0; mv < 4; ++mv )
               {
                  // logical index in the macrocell
                  indexing::Index idx = idx1_cell + vertexdof::logicalIndexOffsetFromVertex( microCell[mv] );
                  // coordinate offset from center
                  coordinateOffset[c][mc][mv] = vertexdof::macrocell::coordinateFromIndex( lvl, *cell, idx ) - center;
                  // micro-vertex on macro-edge, -face, or -cell?
                  auto                  faces_mc_k = vertexdof::macrocell::isOnCellFace( idx, lvl );
                  std::vector< uint_t > intersectingFaces;
                  std::set_intersection( edgeFaces.begin(),
                                         edgeFaces.end(),
                                         faces_mc_k.begin(),
                                         faces_mc_k.end(),
                                         std::back_inserter( intersectingFaces ) );

                  // compute stencil index for edge stencil
                  if ( intersectingFaces.size() >= 2 ) // edge
                  {
                     const auto edgeIdx = indexing::basisConversion( idx, { 0, 1, 2, 3 }, indexingBasis, n );
                     const int  offset  = int( edgeIdx.x() - 1 );
                     if ( offset == 0 )
                     {
                        stencilIndex[c][mc][mv] = vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_C );
                     }
                     else if ( offset == 1 )
                     {
                        stencilIndex[c][mc][mv] = vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_E );
                     }
                     else if ( offset == -1 )
                     {
                        stencilIndex[c][mc][mv] = vertexdof::macroedge::stencilIndexOnEdge( stencilDirection::VERTEX_W );
                     }
                  }
                  else if ( intersectingFaces.size() == 1 ) // face
                  {
                     const auto faceId        = cell->neighborFaces()[intersectingFaces[0]];
                     const auto faceIdxOnEdge = edge->face_index( faceId );
                     // To get the correct indexing basis, we check which one results in a zero entry in the z coordinate.
                     const auto                    firstTestIndexingBasis  = indexingBasis;
                     const std::array< uint_t, 4 > secondTestIndexingBasis = {
                         indexingBasis[0], indexingBasis[1], indexingBasis[3], indexingBasis[2] };
                     const auto edgeIdxFirst  = indexing::basisConversion( idx, { 0, 1, 2, 3 }, firstTestIndexingBasis, n );
                     const auto edgeIdxSecond = indexing::basisConversion( idx, { 0, 1, 2, 3 }, secondTestIndexingBasis, n );
                     const auto faceIdx       = edgeIdxFirst.z() == 0 ? edgeIdxFirst : edgeIdxSecond;
                     const int  offset        = int( faceIdx.x() - 1 );
                     if ( offset == 0 )
                     {
                        stencilIndex[c][mc][mv] =
                            vertexdof::macroedge::stencilIndexOnNeighborFace( stencilDirection::VERTEX_E, faceIdxOnEdge );
                     }
                     else if ( offset == 1 )
                     {
                        stencilIndex[c][mc][mv] =
                            vertexdof::macroedge::stencilIndexOnNeighborFace( stencilDirection::VERTEX_W, faceIdxOnEdge );
                     }
                  }
                  else if ( intersectingFaces.size() == 0 ) // cell
                  {
                     stencilIndex[c][mc][mv] = vertexdof::macroedge::stencilIndexOnNeighborCell( c, n_nbr_faces );
                  }
               }
            }
         }

         // loop over inner vertices on the macro edge
         for ( uint_t i = 1; i < n - 1; ++i )
         {
            auto x = x0 + real_t( i ) * dx;

            auto& stencil = stencils[i - 1];
            stencil.resize( stencil_size );
            std::fill( stencil.begin(), stencil.end(), real_t( 0 ) );

            // loop over neighboring cells
            for ( uint_t c = 0; c < n_nbr_cells; ++c )
            {
               // loop over micro cells
               for ( uint_t mc = 0; mc < n_microCells[c]; ++mc )
               {
                  // coordinates of the micro cell
                  std::array< Point3D, 4 > coords = coordinateOffset[c][mc];
                  for ( auto& coord : coords )
                  {
                     coord += x;
                  }
                  // compute local stiffness matrix
                  form[c].integrateRow( 0, coords, localStiffnessMatrixRow );
                  // assemble stencil
                  for ( uint_t mv = 0; mv < 4; ++mv )
                  {
                     stencil[stencilIndex[c][mc][mv]] += localStiffnessMatrixRow[mv];
                  }
               }
            }
         }
      }
   }

   void precompute_stencil_face_2d( uint_t lvl )
   {
      const auto n     = levelinfo::num_microvertices_per_edge( lvl );
      const auto h     = real_t( 1.0 / ( real_t( n - 1 ) ) );
      const auto n_dof = levelinfo::num_microvertices_per_face_from_width( n - 2 );

      for ( const auto& [faceId, face] : storage_->getFaces() )
      {
         auto& stencils = stencil_face_2d_[faceId][lvl];
         stencils.resize( n_dof );
         form_.setGeometryMap( face->getGeometryMap() );

         // coordinates
         const Point3D x0 = face->getCoordinates()[0];
         const Point3D dx = h * ( face->getCoordinates()[1] - x0 );
         const Point3D dy = h * ( face->getCoordinates()[2] - x0 );

         // coordinate offsets of stencil nbrs
         p1::stencil::StencilData< 2, Point3D > d;
         d[p1::stencil::W]  = -dx;
         d[p1::stencil::E]  = dx;
         d[p1::stencil::N]  = dy;
         d[p1::stencil::S]  = -dy;
         d[p1::stencil::NW] = d[p1::stencil::N] + d[p1::stencil::W];
         d[p1::stencil::SE] = d[p1::stencil::S] + d[p1::stencil::E];

         p1::stencil::StencilData< 2, Point3D >                   coords;
         const std::array< std::array< p1::stencil::Dir, 2 >, 6 > microElements{ { p1::stencil::W, p1::stencil::S },
                                                                                 { p1::stencil::S, p1::stencil::SE },
                                                                                 { p1::stencil::SE, p1::stencil::E },
                                                                                 { p1::stencil::E, p1::stencil::N },
                                                                                 { p1::stencil::N, p1::stencil::NW },
                                                                                 { p1::stencil::NW, p1::stencil::W } };

         // loop over inner vertices on the macro face
         uint_t dof = 0;
         for ( uint_t j = 1; j < n - 2; ++j )
         {
            const Point3D xj = x0 + real_t( j ) * dy;

            for ( uint_t i = 1; i < n - 1 - j; ++i )
            {
               const Point3D x = xj + real_t( i ) * dx;

               for ( uint_t dir = 1; dir < d.size(); ++dir )
               {
                  coords[dir] = x + d[dir];
               }

               // compute local stiffness matrices and add contributions to stencil
               Matrixr< 1, 3 > matrixRow;
               auto&           stencil = stencils[dof];
               stencil.fill( real_t( 0.0 ) );
               for ( auto& el : microElements )
               {
                  form_.integrateRow( 0, { { x, coords[el[0]], coords[el[1]] } }, matrixRow );
                  stencil[p1::stencil::C] += real_t( matrixRow( 0, 0 ) );
                  stencil[el[0]] += real_t( matrixRow( 0, 1 ) );
                  stencil[el[1]] += real_t( matrixRow( 0, 2 ) );
               }
               ++dof;
            }
         }
      }
   }

   void precompute_stencil_face_3d( uint_t lvl )
   {
      const auto n     = levelinfo::num_microvertices_per_edge( lvl );
      const auto h     = real_t( 1.0 / ( real_t( n - 1 ) ) );
      const auto n_dof = levelinfo::num_microvertices_per_face_from_width( n - 2 );

      const indexing::Index idx1( 1, 1, 0 );

      for ( const auto& [faceId, face] : storage_->getFaces() )
      {
         auto& stencils = stencil_face_3d_[faceId][lvl];
         stencils.resize( n_dof );

         // coordinates
         const Point3D x0 = face->getCoordinates()[0];
         const Point3D dx = h * ( face->getCoordinates()[1] - x0 );
         const Point3D dy = h * ( face->getCoordinates()[2] - x0 );

         // coordinate offsets of stencil nbrs
         p1::stencil::StencilData< 3, Point3D > d;

         // logical offsets
         const auto& logical_offsets = offsets_face_3d_.at( faceId );

         // neighbor cells (either 1 or 2)
         const auto  n_nbr_cells = face->getNumNeighborCells();
         const auto& nbrCells    = face->neighborCells();

         // collect required data for each neighbor cell
         std::array< P1Form, 2 >                                 form{ form_, form_ };
         sdt::array< uint_t, 2 >                                 n_microCells{};
         std::array< std::vector< std::array< uint_t, 4 > >, 2 > stencilIndex{};

         for ( uint_t c = 0; c < n_nbr_cells; ++c )
         {
            const Cell* cell = storage_->getCell( nbrCells[c] );
            form[c].setGeometryMap( cell->getGeometryMap() );

            // convert from face to cell coordinate system
            const uint_t                  localFaceId   = cell->getLocalFaceID( faceId );
            const std::array< uint_t, 4 > indexingBasis = algorithms::getMissingIntegersAscending< 3, 4 >(
                { cell->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 ),
                  cell->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 ),
                  cell->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 ) } );
            const auto idx1_cell = indexing::basisConversion( idx1, indexingBasis, { 0, 1, 2, 3 }, n );
            const auto center    = vertexdof::macrocell::coordinateFromIndex( lvl, *cell, idx1_cell );

            // micro-cells associated with the stencil
            const auto microCells = P1Elements::P1Elements3D::getNeighboringElements( idx1_cell, lvl );
            n_microCells[c]       = microCells.size();
            coordinateOffset[c].resize( microCells.size() );
            stencilIndex[c].resize( microCells.size() );

            // iterate over micro-cells
            for ( uint_t mc = 0; mc < microCells.size(); ++mc )
            {
               const auto& microCell = microCells[mc];
               // iterate over the micro-vertices of the micro-cell
               for ( uint_t mv = 0; mv < 4; ++mv )
               {
                  // logical index in the macrocell
                  const indexing::Index idx = idx1_cell + vertexdof::logicalIndexOffsetFromVertex( microCell[mv] );
                  // logical index on the face
                  const auto idx_face = indexing::basisConversion( idx, { 0, 1, 2, 3 }, indexingBasis, n );

                  // stencil direction
                  for ( uint_t dir = 0; dir < logical_offsets.size(); ++dir )
                  {
                     if ( idx_face == logical_offsets[dir] )
                     {
                        stencilIndex[c][mc][mv] = dir;
                        break;
                     }
                  }

                  // coordinate offset from center
                  d[stencilIndex[c][mc][mv]] = vertexdof::macrocell::coordinateFromIndex( lvl, *cell, idx ) - center;
               }
            }
         }

         // assemble stencil
         p1::stencil::StencilData< 3, Point3D > coords;
         Matrixr< 1, 4 >                        localStiffnessMatrixRow;

         // loop over inner vertices on the macro face
         uint_t dof = 0;
         for ( uint_t j = 1; j < n - 2; ++j )
         {
            const Point3D xj = x0 + real_t( j ) * dy;

            for ( uint_t i = 1; i < n - 1 - j; ++i )
            {
               const Point3D x = xj + real_t( i ) * dx;

               for ( uint_t dir = 1; dir < d.size(); ++dir )
               {
                  coords[dir] = x + d[dir];
               }

               auto& stencil = stencils[dof];
               stencil.fill( real_t( 0.0 ) );

               // loop over neighboring cells
               for ( uint_t c = 0; c < n_nbr_cells; ++c )
               {
                  // loop over micro cells
                  for ( uint_t mc = 0; mc < n_microCells[c]; ++mc )
                  {
                     // coordinates of the micro cell
                     std::array< Point3D, 4 > coords_mc;
                     for ( uint_t mv = 0; mv < 4; ++mv )
                     {
                        coords_mc[mv] = coords[stencilIndex[c][mc][mv]];
                     }
                     // compute local stiffness matrix
                     form[c].integrateRow( 0, coords_mc, localStiffnessMatrixRow );
                     // assemble stencil
                     for ( uint_t mv = 0; mv < 4; ++mv )
                     {
                        stencil[stencilIndex[c][mc][mv]] += localStiffnessMatrixRow[mv];
                     }
                  }
               }
               ++dof;
            }
         }
      }
   }

   void precompute_stencil_cell_3d( uint_t lvl )
   {
      const auto n     = levelinfo::num_microvertices_per_edge( lvl );
      const auto h     = real_t( 1.0 / ( real_t( n - 1 ) ) );
      const auto n_dof = levelinfo::num_microvertices_per_cell_from_width( n - 2 );

      std::array< std::array< p1::stencil::Dir, 4 >, 24 > microElements;
      for ( uint_t mc = 0; mc < 24; ++mc )
      {
         for ( uint_t mv = 0; mv < 4; ++mv )
         {
            microElements[mc][mv] = p1::stencil::conversion( P1Elements::P1Elements3D::allCellsAtInnerVertex[mc][mv] )
         }
      }

      for ( const auto& [cellId, cell] : storage_->getCells() )
      {
         auto& stencils = stencil_cell_3d_[cellId][lvl];
         stencils.resize( n_dof );
         form_.setGeometryMap( cell->getGeometryMap() );

         // coordinates
         const Point3D x0 = cell->getCoordinates()[0];
         const Point3D dx = h * ( cell->getCoordinates()[1] - x0 );
         const Point3D dy = h * ( cell->getCoordinates()[2] - x0 );
         const Point3D dz = h * ( cell->getCoordinates()[3] - x0 );

         // coordinate offsets of stencil nbrs
         p1::stencil::StencilData< 3, Point3D > d;
         d[p1::stencil::W]   = -dx;
         d[p1::stencil::E]   = dx;
         d[p1::stencil::N]   = dy;
         d[p1::stencil::S]   = -dy;
         d[p1::stencil::NW]  = d[p1::stencil::N] + d[p1::stencil::W];
         d[p1::stencil::SE]  = d[p1::stencil::S] + d[p1::stencil::E];
         d[p1::stencil::TC]  = dz;
         d[p1::stencil::TW]  = d[p1::stencil::TC] + d[p1::stencil::W];
         d[p1::stencil::TS]  = d[p1::stencil::TC] + d[p1::stencil::S];
         d[p1::stencil::TSE] = d[p1::stencil::TC] + d[p1::stencil::SE];
         d[p1::stencil::BC]  = -dz;
         d[p1::stencil::BE]  = d[p1::stencil::BC] + d[p1::stencil::E];
         d[p1::stencil::BN]  = d[p1::stencil::BC] + d[p1::stencil::N];
         d[p1::stencil::BNW] = d[p1::stencil::BC] + d[p1::stencil::NW];

         p1::stencil::StencilData< 3, Point3D > coords;

         // loop over inner vertices on the macro cell
         uint_t dof = 0;
         for ( uint_t k = 1; k < n - 3; ++k )
         {
            const Point3D xk = x0 + real_t( k ) * dz;
            for ( uint_t j = 1; j < n - 2 - k; ++j )
            {
               const Point3D xj = xk + real_t( j ) * dy;

               for ( uint_t i = 1; i < n - 1 - j - k; ++i )
               {
                  const Point3D x = xj + real_t( i ) * dx;

                  for ( uint_t dir = 1; dir < d.size(); ++dir )
                  {
                     coords[dir] = x + d[dir];
                  }

                  // compute local stiffness matrices and add contributions to stencil
                  Matrixr< 1, 4 > matrixRow;
                  auto&           stencil = stencils[dof];
                  stencil.fill( real_t( 0.0 ) );
                  for ( auto& el : microElements )
                  {
                     form_.integrateRow( 0, { { x, coords[el[1]], coords[el[2]], coords[el[3]] } }, matrixRow );
                     for ( uint_t mv = 0; mv < 4; ++mv )
                     {
                        stencil[el[mv]] += real_t( matrixRow( 0, mv ) );
                     }
                  }
                  ++dof;
               }
            }
         }
      }
   }

   void compute_surrogates_edge_2d( uint_t lvl )
   {
      auto&        lsq = *lsq_interface_[lvl];
      LSQData< 2 > samples;
      for ( uint_t d = 0; d < 7; ++d )
      {
         samples[d].resize( lsq.rows );
      }

      const auto n = levelinfo::num_microvertices_per_edge( lvl );
      const auto h = real_t( 1.0 / ( real_t( n - 1 ) ) );

      P1Form      form_N( form_ );
      Face*       face_S;
      Face*       face_N;
      PrimitiveID vtxId_N, vtxId_S;

      for ( const auto& [edgeId, edge] : storage_->getEdges() )
      {
         // neighbor face
         auto n_nbr_faces = edge->getNumNeighborFaces();
         face_S           = storage_->getFace( edge->neighborFaces()[0] );
         form_.setGeometryMap( face_S->getGeometryMap() );
         vtxId_S = face_S->get_vertex_opposite_to_edge( edgeId );
         if ( n_nbr_faces == 2 )
         {
            face_N = storage_->getFace( edge->neighborFaces()[1] );
            form_N.setGeometryMap( face_N->getGeometryMap() );
            vtxId_N = face_N->get_vertex_opposite_to_edge( edgeId );
         }

         // coordinates
         auto x0 = edge->getCoordinates()[0];
         auto x1 = edge->getCoordinates()[1];
         auto dx = h * edge->getDirection();
         // coordinate offsets of stencil nbrs
         p1::stencil::StencilData< 2, Point3D > d;
         d[p1::stencil::W] = -dx;
         d[p1::stencil::E] = dx;
         if ( n_nbr_faces == 2 )
         {
            auto coord_N      = storage_->getVertex( vtxId_N )->getCoordinates();
            d[p1::stencil::N] = h * ( coord_N - x0 );
         }
         auto coord_S       = storage_->getVertex( vtxId_S )->getCoordinates();
         d[p1::stencil::S]  = h * ( coord_S - x1 );
         d[p1::stencil::NW] = d[p1::stencil::N] + d[p1::stencil::W];
         d[p1::stencil::SE] = d[p1::stencil::S] + d[p1::stencil::E];

         // loop over inner vertices on the macro edge
         auto it = lsq.samplingIterator();
         while ( it != it.end() )
         {
            // compute local stiffness matrices and add contributions to stencil
            Stencil< 2 > stencil {}
            Matrixr< 1, 3 > matrixRow;
            auto            x = x0 + real_t( it.i() ) * dx;
            form_.integrateRow( 0, { { x, x + d[p1::stencil::W], x + d[p1::stencil::S] } }, matrixRow );
            stencil[p1::stencil::C] = real_t( matrixRow( 0, 0 ) );
            stencil[p1::stencil::W] = real_t( matrixRow( 0, 1 ) );
            stencil[p1::stencil::S] = real_t( matrixRow( 0, 2 ) );
            form_.integrateRow( 0, { { x, x + d[p1::stencil::S], x + d[p1::stencil::SE] } }, matrixRow );
            stencil[p1::stencil::C] += real_t( matrixRow( 0, 0 ) );
            stencil[p1::stencil::S] += real_t( matrixRow( 0, 1 ) );
            stencil[p1::stencil::SE] = real_t( matrixRow( 0, 2 ) );
            form_.integrateRow( 0, { { x, x + d[p1::stencil::SE], x + d[p1::stencil::E] } }, matrixRow );
            stencil[p1::stencil::C] += real_t( matrixRow( 0, 0 ) );
            stencil[p1::stencil::SE] += real_t( matrixRow( 0, 1 ) );
            stencil[p1::stencil::E] = real_t( matrixRow( 0, 2 ) );
            if ( n_nbr_faces == 2 )
            {
               form_N.integrateRow( 0, { { x, x + d[p1::stencil::E], x + d[p1::stencil::N] } }, matrixRow );
               stencil[p1::stencil::C] += real_t( matrixRow( 0, 0 ) );
               stencil[p1::stencil::E] += real_t( matrixRow( 0, 1 ) );
               stencil[p1::stencil::N] = real_t( matrixRow( 0, 2 ) );
               form_N.integrateRow( 0, { { x, x + d[p1::stencil::N], x + d[p1::stencil::NW] } }, matrixRow );
               stencil[p1::stencil::C] += real_t( matrixRow( 0, 0 ) );
               stencil[p1::stencil::N] += real_t( matrixRow( 0, 1 ) );
               stencil[p1::stencil::NW] = real_t( matrixRow( 0, 2 ) );
               form_N.integrateRow( 0, { { x, x + d[p1::stencil::NW], x + d[p1::stencil::W] } }, matrixRow );
               stencil[p1::stencil::C] += real_t( matrixRow( 0, 0 ) );
               stencil[p1::stencil::NW] += real_t( matrixRow( 0, 1 ) );
               stencil[p1::stencil::W] += real_t( matrixRow( 0, 2 ) );
            }

            // add data sample
            for ( uint_t d = 0; d < stencil.size(); ++d )
            {
               samples[d][it()] = stencil[d];
            }
            ++it;
         }
         for ( uint_t i = 1; i < n - 1; ++i ) {}

         // fit stencil polynomials
         PolyStencil< 2, 1 >& surrogate = surrogate_edge_2d_[edgeId][lvl];
         for ( uint_t d = 0; d < stencil.size(); ++d )
         {
            lsq.setRHS( samples[d] );
            auto& coeffs      = lsq.solve();
            surrogate.set_coefficients(d, coeffs);
         }
      }
   }

   void compute_surrogates_face_2d( uint_t lvl )
   {
      // todo
   }

   void compute_surrogates_face_3d( uint_t lvl )
   {
      // todo
   }

   void compute_surrogates_cell_3d( uint_t lvl )
   {
      // todo
   }

   void compute_offsets_for_face_stencil_3d()
   {
      /* only the first 7 directions are the same for all faces. The directions in the cell interior depend on the
         local orientation of the face w.r.t. the cell. */

      const uint_t          lvl = 2; // make sure there are interior vertices
      const auto            n   = levelinfo::num_microvertices_per_edge( lvl );
      const indexing::Index idx1( 1, 1, 0 );

      for ( const auto& [faceId, face] : storage_->getFaces() )
      {
         p1::stencil::StencilData< 3, indexing::Index > logical_offsets = p1::stencil::offset;

         // neighbor cells (either 1 or 2)
         const auto  n_nbr_cells = face->getNumNeighborCells();
         const auto& nbrCells    = face->neighborCells();

         for ( uint_t c = 0; c < n_nbr_cells; ++c )
         {
            const Cell* cell = storage_->getCell( nbrCells[c] );

            // convert from face to cell coordinate system
            const uint_t                  localFaceId   = cell->getLocalFaceID( faceId );
            const std::array< uint_t, 4 > indexingBasis = algorithms::getMissingIntegersAscending< 3, 4 >(
                { cell->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 0 ),
                  cell->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 1 ),
                  cell->getFaceLocalVertexToCellLocalVertexMaps().at( localFaceID ).at( 2 ) } );
            const auto idx1_cell = indexing::basisConversion( idx1, indexingBasis, { 0, 1, 2, 3 }, n );

            // micro-cells associated with the stencil
            const auto microCells = P1Elements::P1Elements3D::getNeighboringElements( idx1_cell, lvl );

            // iterate over micro-cells
            uint_t n_offsets_assigned = 7 + 4 * c;
            for ( const auto& microCell : microCells )
            {
               // iterate over the micro-vertices of the micro-cell
               for ( auto& mv : microCell )
               {
                  // logical index in the macrocell
                  const indexing::Index idx = idx1_cell + vertexdof::logicalIndexOffsetFromVertex( mv );
                  // logical index on the face
                  const auto idx_face = indexing::basisConversion( idx, { 0, 1, 2, 3 }, indexingBasis, n );

                  // stencil direction
                  if ( idx_face.z() != 0 ) // micro-vertex in interior of cell
                  {
                     // check whether idx_face is already associated with some direction
                     bool not_registered_yet =
                         std::find( logical_offsets.begin(), logical_offsets.begin() + n_offsets_assigned, idx_face ) ==
                         logical_offsets.begin() + n_offsets_assigned;

                     // define association between stencil::Dir and idx_face
                     if ( not_registered_yet )
                     {
                        logical_offsets[n_offsets_assigned] = idx_face;
                        ++n_offsets_assigned;
                     }
                  }
               }
            }
         }
         offsets_face_3d_[faceId] = logical_offsets;
      }
   }

   void apply_edge_precomputed_2d( std::shared_ptr< hyteg::Edge > edge, uint_t lvl, const real_t* srcData, real_t* dstData )
   {
      const auto n           = levelinfo::num_microvertices_per_edge( lvl );
      const auto n_nbr_faces = edge->getNumNeighborFaces();
      const auto stencilSize = 3 + 2 * n_nbr_faces;

      // indices of neighboring DoF
      DofIdx dofIdx{};
      for ( int d = 0; d < stencilSize; ++d )
      {
         dofIdx[d] = vertexdof::macroedge::indexFromVertex( lvl, 1, p1::stencil::backConversion[d] );
      }

      const auto& stencils = stencil_edge_2d_.at( edge->getID() )[lvl];

      // loop over inner vertices on the macro edge
      for ( uint_t i = 1; i < n - 1; ++i )
      {
         const auto& stencil = stencils[i - 1];
         const auto  dstIdx  = dofIdx[p1::stencil::C];

         if ( updateType == Replace )
         {
            dstData[i] = real_c( 0 );
         } // else updateType == Add

         // apply stencil
         for ( int d = 0; d < stencilSize; ++d )
         {
            dstData[dstIdx] += stencil[d] * srcData[dofIdx[d]];
            ++dofIdx[d];
         }
      }
   }

   void apply_edge_precomputed_3d( std::shared_ptr< hyteg::Edge > edge, uint_t lvl, const real_t* srcData, real_t* dstData )
   {
      const auto n           = levelinfo::num_microvertices_per_edge( lvl );
      const auto n_nbr_faces = edge->getNumNeighborFaces();
      const auto n_nbr_cells = edge->getNumNeighborCells();
      const auto stencilSize = vertexDoFMacroEdgeStencilMemorySize( lvl, edge );

      constexpr auto C = stencilDirection::VERTEX_C;
      constexpr auto W = stencilDirection::VERTEX_W;
      constexpr auto E = stencilDirection::VERTEX_E;

      // indices of neighboring DoF
      std::vector< walberla::uint_t > dofIdx( stencilSize );
      // for the 3D edge the center is not associated with stencil[0], since, for simplicity, we keep the old ordering here.
      uint_t centerIdx = vertexdof::macroedge::stencilIndexOnEdge( C );
      for ( auto& sd : { C, W, E } )
      {
         dofIdx[vertexdof::macroedge::stencilIndexOnEdge( sd )] = //
             vertexdof::macroedge::indexFromVertex( lvl, 1, sd );
      }
      for ( uint_t faceIdx = 0; faceIdx < n_nbr_faces; ++faceIdx )
      {
         for ( auto& sd : { W, E } )
         {
            dofIdx[vertexdof::macroedge::stencilIndexOnNeighborFace( sd, faceIdx )] =
                vertexdof::macroedge::indexFromVertexOnNeighborFace( lvl, 1, faceIdx, sd );
         }
      }
      for ( uint_t cellIdx = 0; cellIdx < n_nbr_cells; ++cellIdx )
      {
         dofIdx[vertexdof::macroedge::stencilIndexOnNeighborCell( cellIdx, n_nbr_faces )] =
             vertexdof::macroedge::indexFromVertexOnNeighborCell( lvl, 1, cellIdx, n_nbr_faces );
      }

      const auto& stencils = stencil_edge_3d_.at( edge->getID() )[lvl];

      // loop over inner vertices on the macro edge
      for ( uint_t i = 1; i < n - 1; ++i )
      {
         const auto& stencil = stencils[i - 1];
         const auto  dstIdx  = dofIdx[centerIdx];

         if ( updateType == Replace )
         {
            dstData[i] = real_c( 0 );
         } // else updateType == Add

         // apply stencil
         for ( int d = 0; d < stencilSize; ++d )
         {
            dstData[dstIdx] += stencil[d] * srcData[dofIdx[d]];
            ++dofIdx[d];
         }
      }
   }

   void apply_face_precomputed_2d( std::shared_ptr< hyteg::Face > face, uint_t lvl, const real_t* srcData, real_t* dstData )
   {
      const auto     n           = levelinfo::num_microvertices_per_edge( lvl );
      constexpr auto stencilSize = 7;

      const auto& stencils = stencil_face_2d_.at( face->getID() )[lvl];

      // loop over inner vertices on the macro face
      uint_t dof = 0;
      for ( uint_t j = 1; j < n - 2; ++j )
      {
         // indices of neighboring DoF
         DofIdx dofIdx{};
         for ( int d = 0; d < stencilSize; ++d )
         {
            dofIdx[d] = vertexdof::macroface::indexFromVertex( lvl, 1, j, p1::stencil::backConversion[d] );
         }

         for ( uint_t i = 1; i < n - 1 - j; ++i )
         {
            const auto& stencil = stencils[dof];
            const auto  dstIdx  = dofIdx[p1::stencil::C];

            if ( updateType == Replace )
            {
               dstData[i] = real_c( 0 );
            } // else updateType == Add

            // apply stencil
            for ( int d = 0; d < stencilSize; ++d )
            {
               dstData[dstIdx] += stencil[d] * srcData[dofIdx[d]];
               ++dofIdx[d];
            }

            ++dof;
         }
      }
   }

   void apply_face_precomputed_3d( std::shared_ptr< hyteg::Face > face, uint_t lvl, const real_t* srcData, real_t* dstData )
   {
      const auto  n           = levelinfo::num_microvertices_per_edge( lvl );
      const auto  n_nbr_cells = face->getNumNeighborCells();
      const auto  stencilSize = 7 + 4 * n_nbr_cells;
      const auto& offsets     = offsets_face_3d_.at( faceId );

      const auto& stencils = stencil_face_3d_.at( face->getID() )[lvl];

      // loop over inner vertices on the macro face
      uint_t dof = 0;
      for ( uint_t j = 1; j < n - 2; ++j )
      {
         // indices of neighboring DoF
         DofIdx dofIdx{};
         for ( int d = 0; d < 7; ++d )
         {
            // neighbors on face
            dofIdx[d] = vertexdof::macroface::index( lvl, 1 + offsets[d].x(), j + offsets[d].y() );
         }
         for ( int d = 7; d < 11; ++d )
         {
            // neighbors on first nbr cell
            dofIdx[d] = vertexdof::macroface::index( lvl, 1 + offsets[d].x(), j + offsets[d].y(), 0 );
         }
         for ( int d = 11; d < stencilSize; ++d )
         {
            // neighbors on second nbr cell
            dofIdx[d] = vertexdof::macroface::index( lvl, 1 + offsets[d].x(), j + offsets[d].y(), 1 );
         }

         for ( uint_t i = 1; i < n - 1 - j; ++i )
         {
            const auto& stencil = stencils[dof];
            const auto  dstIdx  = dofIdx[p1::stencil::C];

            if ( updateType == Replace )
            {
               dstData[i] = real_c( 0 );
            } // else updateType == Add

            // apply stencil
            for ( int d = 0; d < stencilSize; ++d )
            {
               dstData[dstIdx] += stencil[d] * srcData[dofIdx[d]];
               ++dofIdx[d];
            }

            ++dof;
         }
      }
   }

   void apply_cell_precomputed_3d( std::shared_ptr< hyteg::Cell > cell, uint_t lvl, const real_t* srcData, real_t* dstData )
   {
      const auto     n           = levelinfo::num_microvertices_per_edge( lvl );
      constexpr auto stencilSize = 15;

      const auto& stencils = stencil_cell_3d_.at( cell->getID() )[lvl];

      // loop over inner vertices on the macro cell
      uint_t dof = 0;
      for ( uint_t k = 1; k < n - 3; ++k )
      {
         for ( uint_t j = 1; j < n - 2 - k; ++j )
         {
            // indices of neighboring DoF
            DofIdx dofIdx{};
            for ( int d = 0; d < stencilSize; ++d )
            {
               dofIdx[d] = vertexdof::macrocell::indexFromVertex( lvl, 1, j, k, p1::stencil::backConversion[d] );
            }

            for ( uint_t i = 1; i < n - 1 - j - k; ++i )
            {
               const auto& stencil = stencils[dof];
               const auto  dstIdx  = dofIdx[p1::stencil::C];

               if ( updateType == Replace )
               {
                  dstData[i] = real_c( 0 );
               } // else updateType == Add

               // apply stencil
               for ( int d = 0; d < stencilSize; ++d )
               {
                  dstData[dstIdx] += stencil[d] * srcData[dofIdx[d]];
                  ++dofIdx[d];
               }

               ++dof;
            }
         }
      }
   }

   void apply_edge_surrogate_2d( std::shared_ptr< hyteg::Edge > edge, uint_t lvl, const real_t* srcData, real_t* dstData )
   {
      const auto n           = levelinfo::num_microvertices_per_edge( lvl );
      const auto n_nbr_faces = edge->getNumNeighborFaces();
      const auto stencilSize = 3 + 2 * n_nbr_faces;

      // indices of neighboring DoF
      DofIdx dofIdx{};
      for ( int d = 0; d < stencilSize; ++d )
      {
         dofIdx[d] = vertexdof::macroedge::indexFromVertex( lvl, 1, p1::stencil::backConversion[d] );
      }

      const PolyStencil< 2, 1 >& surrogate = surrogate_edge_2d_.at( edge->getID() )[lvl];
      const PolyDomain           X( lvl );

      // loop over inner vertices on the macro edge
      for ( uint_t i = 1; i < n - 1; ++i )
      {
         // evaluate polynomial
         const auto x = X[i];
         surrogate.eval( x );
         const auto& stencil = surrogate.px();

         const auto dstIdx = dofIdx[p1::stencil::C];

         if ( updateType == Replace )
         {
            dstData[i] = real_c( 0 );
         } // else updateType == Add

         // apply stencil
         for ( int d = 0; d < stencilSize; ++d )
         {
            dstData[dstIdx] += stencil[d] * srcData[dofIdx[d]];
            ++dofIdx[d];
         }
      }
   }

   void apply_face_surrogate_2d( std::shared_ptr< hyteg::Face > face, uint_t lvl, const real_t* srcData, real_t* dstData )
   {
      const auto     n           = levelinfo::num_microvertices_per_edge( lvl );
      constexpr auto stencilSize = 7;

      const PolyStencil< 2, 2 >& surrogate = surrogate_face_2d_.at( face->getID() )[lvl];
      const PolyDomain           X( lvl );

      // loop over inner vertices on the macro face
      for ( uint_t j = 1; j < n - 2; ++j )
      {
         // indices of neighboring DoF
         DofIdx dofIdx{};
         for ( int d = 0; d < stencilSize; ++d )
         {
            dofIdx[d] = vertexdof::macroface::indexFromVertex( lvl, 1, j, p1::stencil::backConversion[d] );
         }

         // restrict polynomial to 1D
         const auto y = X[j];
         surrogate.fix_y( y );

         for ( uint_t i = 1; i < n - 1 - j; ++i )
         {
            // evaluate polynomial
            const auto x = X[i];
            surrogate.eval( x );
            const auto& stencil = surrogate.px();

            const auto dstIdx = dofIdx[p1::stencil::C];

            if ( updateType == Replace )
            {
               dstData[i] = real_c( 0 );
            } // else updateType == Add

            // apply stencil
            for ( int d = 0; d < stencilSize; ++d )
            {
               dstData[dstIdx] += stencil[d] * srcData[dofIdx[d]];
               ++dofIdx[d];
            }
         }
      }
   }

   void apply_face_surrogate_3d( std::shared_ptr< hyteg::Face > face, uint_t lvl, const real_t* srcData, real_t* dstData )
   {
      const auto  n           = levelinfo::num_microvertices_per_edge( lvl );
      const auto  n_nbr_cells = face->getNumNeighborCells();
      const auto  stencilSize = 7 + 4 * n_nbr_cells;
      const auto& offsets     = offsets_face_3d_.at( faceId );

      const PolyStencil< 3, 2 >& surrogate = surrogate_face_3d_.at( face->getID() )[lvl];
      const PolyDomain           X( lvl );

      // loop over inner vertices on the macro face
      for ( uint_t j = 1; j < n - 2; ++j )
      {
         // indices of neighboring DoF
         DofIdx dofIdx{};
         for ( int d = 0; d < 7; ++d )
         {
            // neighbors on face
            dofIdx[d] = vertexdof::macroface::index( lvl, 1 + offsets[d].x(), j + offsets[d].y() );
         }
         for ( int d = 7; d < 11; ++d )
         {
            // neighbors on first nbr cell
            dofIdx[d] = vertexdof::macroface::index( lvl, 1 + offsets[d].x(), j + offsets[d].y(), 0 );
         }
         for ( int d = 11; d < stencilSize; ++d )
         {
            // neighbors on second nbr cell
            dofIdx[d] = vertexdof::macroface::index( lvl, 1 + offsets[d].x(), j + offsets[d].y(), 1 );
         }

         // restrict polynomial to 1D
         const auto y = X[j];
         surrogate.fix_y( y );

         for ( uint_t i = 1; i < n - 1 - j; ++i )
         {
            // evaluate polynomial
            const auto x = X[i];
            surrogate.eval( x );
            const auto& stencil = surrogate.px();

            const auto dstIdx = dofIdx[p1::stencil::C];

            if ( updateType == Replace )
            {
               dstData[i] = real_c( 0 );
            } // else updateType == Add

            // apply stencil
            for ( int d = 0; d < stencilSize; ++d )
            {
               dstData[dstIdx] += stencil[d] * srcData[dofIdx[d]];
               ++dofIdx[d];
            }
         }
      }
   }

   void apply_cell_surrogate_3d( std::shared_ptr< hyteg::Cell > cell, uint_t lvl, const real_t* srcData, real_t* dstData )
   {
      const auto     n           = levelinfo::num_microvertices_per_edge( lvl );
      constexpr auto stencilSize = 15;

      const PolyStencil< 3, 3 >& surrogate = surrogate_cell_3d_.at( cell->getID() )[lvl];
      const PolyDomain           X( lvl );

      // loop over inner vertices on the macro cell
      for ( uint_t k = 1; k < n - 3; ++k )
      {
         // restrict polynomial to 2D
         const auto z = X[k];
         surrogate.fix_z( z );

         for ( uint_t j = 1; j < n - 2 - k; ++j )
         {
            // indices of neighboring DoF
            DofIdx dofIdx{};
            for ( int d = 0; d < stencilSize; ++d )
            {
               dofIdx[d] = vertexdof::macrocell::indexFromVertex( lvl, 1, j, k, p1::stencil::backConversion[d] );
            }

            // restrict polynomial to 1D
            const auto y = X[j];
            surrogate.fix_y( y );

            for ( uint_t i = 1; i < n - 1 - j - k; ++i )
            {
               // evaluate polynomial
               const auto x = X[i];
               surrogate.eval( x );
               const auto& stencil = surrogate.px();

               const auto dstIdx = dofIdx[p1::stencil::C];

               if ( updateType == Replace )
               {
                  dstData[i] = real_c( 0 );
               } // else updateType == Add

               // apply stencil
               for ( int d = 0; d < stencilSize; ++d )
               {
                  dstData[dstIdx] += stencil[d] * srcData[dofIdx[d]];
                  ++dofIdx[d];
               }
            }
         }
      }
   }

   //! todo remove everything from here ========>

   /* interpolate polynomials
   */
   void interpolate2D( uint_t polyDegree, uint_t maxInterpolationLevel )
   {
      for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
      {
         const uint_t interpolationLevel = std::min( level, maxInterpolationLevel );
         const uint_t lvlDiff            = level - interpolationLevel;
         const uint_t rowsizeY           = levelinfo::num_microvertices_per_edge( interpolationLevel );
         const real_t h_il               = 1.0 / ( walberla::real_c( rowsizeY - 1 ) );

         for ( auto& it : storage_->getFaces() )
         {
            Face& face          = *it.second;
            auto  stencilMemory = face.getData( faceStencilID_ )->getPointer( level );
            auto& stencilPoly   = face.getData( facePolyID_ )->getData( level );

            Interpolator_face interpolator;

            assemble_variableStencil_face_init( face, level );

            // initialize polynomials
            // auto stencilSize   = face.getData(faceStencilID_)->getSize(level); // always returns 27!
            auto stencilSize = faceStencilSize2D;

            for ( uint_t c = 0; c < stencilSize; ++c )
            {
               stencilPoly.push_back( Poly2D( polyDegree ) );
               interpolator.push_back( Interpolator2D( polyDegree, interpolationLevel, level < maxLevel_ ) );
            }

            // add sample points
            uint_t  rowsizeX;
            Point2D x;

            for ( uint_t j = 1; j < rowsizeY - 2; ++j )
            {
               x[1]     = j * h_il;
               rowsizeX = rowsizeY - j;

               for ( uint_t i = 1; i < rowsizeX - 1; ++i )
               {
                  x[0] = i * h_il;

                  assemble_variableStencil_face( stencilMemory, i << lvlDiff, j << lvlDiff );

                  for ( uint_t c = 0; c < stencilSize; ++c )
                  {
                     interpolator[c].addInterpolationPoint( x, stencilMemory[c] );
                  }
               }
            }

            // find polynomials by L2 fit
            for ( uint_t c = 0; c < stencilSize; ++c )
            {
               interpolator[c].interpolate( stencilPoly[c] );
            }
         }
      }

      // initialize polynomial evaluator
      // for ( auto& it : storage_->getFaces() )
      // {
      // auto stencilSize   = it.second->getData(faceStencilID_)->getSize(maxLevel_); // always returns 27!
      auto stencilSize = faceStencilSize2D;

      for ( uint_t c = 0; c < stencilSize; ++c )
      {
         facePolyEvaluator_.push_back( Polynomial2DEvaluator( polyDegree ) );
      }

      //    break; // we use the same evaluator for all faces
      // }
   }

   /* interpolate polynomials
   */
   void interpolate3D( uint_t polyDegree, uint_t maxInterpolationLevel )
   {
      for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
      {
         // skip level 0 (no interior points)
         if ( level == 0 )
            continue;

         const uint_t interpolationLevel = std::min( level, maxInterpolationLevel );
         const uint_t lvlDiff            = level - interpolationLevel;
         const uint_t rowsizeZ           = levelinfo::num_microvertices_per_edge( interpolationLevel );
         const real_t h_il               = 1.0 / ( walberla::real_c( rowsizeZ - 1 ) );

         for ( const auto& it : storage_->getCells() )
         {
            Cell& cell          = *it.second;
            auto& stencilMemory = cell.getData( cellStencilID_ )->getData( level );
            auto& stencilPoly   = cell.getData( cellPolyID_ )->getData( level );

            Interpolator_cell interpolator;

            assemble_variableStencil_cell_init( cell, level );

            // initialize polynomials
            assemble_variableStencil_cell( stencilMemory, 1, 1, 1 );

            for ( auto& [idx, val] : stencilMemory )
            {
               stencilPoly.insert_or_assign( idx, Poly3D( polyDegree ) );
               interpolator.insert_or_assign( idx, Interpolator3D( polyDegree, interpolationLevel, level < maxLevel_ ) );
            }

            // add sample points
            uint_t  rowsizeY, rowsizeX;
            Point3D x;

            for ( uint_t k = 1; k < rowsizeZ - 3; ++k )
            {
               x[2]     = k * h_il;
               rowsizeY = rowsizeZ - k;

               for ( uint_t j = 1; j < rowsizeY - 2; ++j )
               {
                  x[1]     = j * h_il;
                  rowsizeX = rowsizeY - j;

                  for ( uint_t i = 1; i < rowsizeX - 1; ++i )
                  {
                     x[0] = i * h_il;

                     assemble_variableStencil_cell( stencilMemory, i << lvlDiff, j << lvlDiff, k << lvlDiff );

                     for ( auto& [idx, val] : stencilMemory )
                     {
                        interpolator[idx].addInterpolationPoint( x, val );
                     }
                  }
               }
            }

            // find polynomials by L2 fit
            for ( auto& [idx, interp] : interpolator )
            {
               interp.interpolate( stencilPoly[idx] );
            }
         }
      }

      // initialize polynomial evaluator
      for ( auto& it : storage_->getCells() )
      {
         auto& stencilMemory = it.second->getData( cellStencilID_ )->getData( maxLevel_ );

         for ( auto& [idx, val] : stencilMemory )
         {
            cellPolyEvaluator_.insert_or_assign( idx, Polynomial3DEvaluator( polyDegree ) );
         }

         break; // we use the same evaluator for all cells
      }
   }

   /// stencil assembly ///////////

   /* Initialize assembly of variable edge stencil.
      Will be called before iterating over edge whenever the stencil is applied.
   */
   inline void assemble_stencil_edge_init( Edge& edge, const uint_t level ) const
   {
      assemble_variableStencil_edge_init( edge, level );
   }

   /* Assembly of edge stencil.
      Will be called before stencil is applied to a particular edge-DoF.
   */
   inline void assemble_stencil_edge( real_t* edge_stencil, const uint_t i ) const
   {
      assemble_variableStencil_edge( edge_stencil, i );
   }

   /* Initialize assembly of face stencil.
      Will be called before iterating over face whenever the stencil is applied.
   */
   inline void assemble_stencil_face_init( Face& face, const uint_t level ) const
   {
      if ( storage_->hasGlobalCells() )
      {
         assemble_variableStencil_face_init( face, level );
      }
      else
      {
         h_ = 1.0 / ( walberla::real_c( levelinfo::num_microvertices_per_edge( level ) - 1 ) );

         auto& stencilPoly = face.getData( facePolyID_ )->getData( level );

         for ( uint_t c = 0; c < facePolyEvaluator_.size(); ++c )
         {
            facePolyEvaluator_[c].setPolynomial( stencilPoly[c] );
         }
      }
   }

   inline void assemble_stencil_face_init_y( const uint_t j ) const
   {
      if ( !( storage_->hasGlobalCells() ) )
      {
         real_t y = h_ * j;

         for ( auto& evaluator : facePolyEvaluator_ )
         {
            evaluator.setY( y );
            if constexpr ( USE_INCREMENTAL_EVAL )
            {
               evaluator.setStartX( 0, h_ );
            }
         }
      }
   }

   /* Assembly of face stencil.
      Will be called before stencil is applied to a particular face-DoF of a 2d domain.
   */
   inline void assemble_stencil_face( real_t* face_stencil, const uint_t i, const uint_t j ) const
   {
      real_t x = h_ * i;

      for ( uint_t c = 0; c < facePolyEvaluator_.size(); ++c )
      {
         if constexpr ( USE_INCREMENTAL_EVAL )
         {
            face_stencil[c] = facePolyEvaluator_[c].incrementEval();
         }
         else
         {
            face_stencil[c] = facePolyEvaluator_[c].evalX( x );
         }
      }
   }

   /* Assembly of face stencil.
      Will be called before stencil is applied to a particular face-DoF of a 3D domain.
   */
   inline void assemble_stencil_face3D( vertexdof::macroface::StencilMap_T& face_stencil, const uint_t i, const uint_t j ) const
   {
      assemble_variableStencil_face3D( face_stencil, i, j );
   }

   /* Initialize assembly of cell stencil.
      Will be called before iterating over cell whenever the stencil is applied.
   */
   inline void assemble_stencil_cell_init( Cell& cell, const uint_t level ) const
   {
      h_ = 1.0 / ( walberla::real_c( levelinfo::num_microvertices_per_edge( level ) - 1 ) );

      auto& stencilPoly = cell.getData( cellPolyID_ )->getData( level );

      for ( auto& [idx, evaluator] : cellPolyEvaluator_ )
      {
         evaluator.setPolynomial( stencilPoly[idx] );
      }
   }

   inline void assemble_stencil_cell_init_z( const uint_t k ) const
   {
      real_t z = h_ * k;

      for ( auto& [idx, evaluator] : cellPolyEvaluator_ )
      {
         evaluator.setZ( z );
      }
   }

   inline void assemble_stencil_cell_init_y( const uint_t j ) const
   {
      real_t y = h_ * j;

      for ( auto& [idx, evaluator] : cellPolyEvaluator_ )
      {
         evaluator.setY( y );
         if constexpr ( USE_INCREMENTAL_EVAL )
         {
            evaluator.setStartX( 0, h_ );
         }
      }
   }

   /* Assembly of cell stencil.
      Will be called before stencil is applied to a particular cell-DoF.
   */
   inline void assemble_stencil_cell( vertexdof::macrocell::StencilMap_T& cell_stencil,
                                      const uint_t                        i,
                                      const uint_t                        j,
                                      const uint_t                        k ) const
   {
      real_t x = h_ * i;

      for ( auto& [idx, evaluator] : cellPolyEvaluator_ )
      {
         if constexpr ( USE_INCREMENTAL_EVAL )
         {
            cell_stencil[idx] = evaluator.incrementEval();
         }
         else
         {
            cell_stencil[idx] = evaluator.evalX( x );
         }
      }
   }

   //! ========> until here

   void computeDiagonalOperatorValues( bool invert )
   {
      auto& diagonal = ( invert ) ? inverseDiagonalValues_ : diagonalValues_;
      if ( !diagonal )
      {
         diagonal = std::make_shared< P1Function< ValueType > >(
             ( invert ) ? "inverse diagonal entries" : "diagonal entries", storage_, minLevel_, maxLevel_ );
      }

      const int dim = ( storage_->hasGlobalCells() ) ? 3 : 2;

      for ( uint_t lvl = minLevel_; lvl <= maxLevel_; ++lvl )
      {
         const auto n = levelinfo::num_microvertices_per_edge( lvl );

         for ( const auto& [vtxId, vertex] : storage_->getVertices() )
         {
            auto        diagData = vertex->getData( diagonal->getVertexDataID() )->getPointer( lvl );
            const auto& stencil  = stencil_vtx_.at( vtxId )[lvl][0];
            diagData[0]          = stencil[0];
         }

         if ( lvl >= 1 )
         {
            for ( const auto& [edgeId, edge] : storage_->getEdges() )
            {
               auto diagData = edge->getData( diagonal->getEdgeDataID() )->getPointer( lvl );

               if ( dim == 2 )
               {
                  if ( lvl < min_lvl_for_surrogate )
                  {
                     assemble_diagonalOperator_edge_precomputed_2d( edge, lvl, diagData );
                  }
                  else
                  {
                     assemble_diagonalOperator_edge_surrogate_2d( edge, lvl, diagData );
                  }
               }
               else // dim == 3
               {
                  assemble_diagonalOperator_edge_precomputed_3d( edge, lvl, diagData );
               }
            }
         }

         if ( lvl >= 2 )
         {
            for ( const auto& [faceId, face] : storage_->getFaces() )
            {
               auto diagData = face->getData( diagonal->getFaceDataID() )->getPointer( lvl );

               if ( dim == 2 )
               {
                  if ( lvl < min_lvl_for_surrogate )
                  {
                     assemble_diagonalOperator_face_precomputed_2d( face, lvl, diagData );
                  }
                  else
                  {
                     assemble_diagonalOperator_face_surrogate_2d( face, lvl, diagData );
                  }
               }
               else // dim == 3
               {
                  if ( lvl < min_lvl_for_surrogate )
                  {
                     assemble_diagonalOperator_face_precomputed_3d( face, lvl, diagData );
                  }
                  else
                  {
                     assemble_diagonalOperator_face_surrogate_3d( face, lvl, diagData );
                  }
               }
            }
         }

         if ( lvl >= 2 )
         {
            for ( const auto& [cellId, cell] : storage_->getCells() )
            {
               if ( lvl < min_lvl_for_surrogate )
               {
                  assemble_diagonalOperator_cell_precomputed_3d( cell, lvl, diagData );
               }
               else
               {
                  assemble_diagonalOperator_cell_surrogate_3d( cell, lvl, diagData );
               }
            }
         }

         if ( invert )
         {
            diagonal->invertElementwise( lvl, All, false );
         }
      }
   }

   void assemble_diagonalOperator_edge_precomputed_2d( std::shared_ptr< hyteg::Edge > edge, uint_t lvl, const real_t* diagData )
   {
      auto        diagIdx  = vertexdof::macroedge::index( lvl, 1 );
      const auto& stencils = stencil_edge_2d_.at( edge->getID() )[lvl];
      for ( uint_t i = 1; i < n - 1; ++i )
      {
         diagData[diagIdx] = stencils[i - 1][p1::stencil::C];
         ++diagIdx;
      }
   }

   void assemble_diagonalOperator_edge_precomputed_3d( std::shared_ptr< hyteg::Edge > edge, uint_t lvl, const real_t* diagData )
   {
      auto        diagIdx  = vertexdof::macroedge::index( lvl, 1 );
      const auto  c        = vertexdof::macroedge::stencilIndexOnEdge( hyteg::stencilDirection::VERTEX_C );
      const auto& stencils = stencil_edge_3d_.at( edgeId )[lvl];
      for ( uint_t i = 1; i < n - 1; ++i )
      {
         diagData[diagIdx] = stencils[i - 1][c];
         ++diagIdx;
      }
   }

   void assemble_diagonalOperator_face_precomputed_2d( std::shared_ptr< hyteg::Face > face, uint_t lvl, const real_t* diagData )
   {
      const auto& stencils = stencil_face_2d_.at( face->getID() )[lvl];

      uint_t dof = 0;
      for ( uint_t j = 1; j < n - 2; ++j )
      {
         auto diagIdx = vertexdof::macroface::index( lvl, 1, j );

         for ( uint_t i = 1; i < n - 1 - j; ++i )
         {
            diagData[diagIdx] = stencils[dof][p1::stencil::C];
            ++diagIdx;
            ++dof;
         }
      }
   }

   void assemble_diagonalOperator_face_precomputed_3d( std::shared_ptr< hyteg::Face > face, uint_t lvl, const real_t* diagData )
   {
      const auto& stencils = stencil_face_3d_.at( face->getID() )[lvl];

      uint_t dof = 0;
      for ( uint_t j = 1; j < n - 2; ++j )
      {
         auto diagIdx = vertexdof::macroface::index( lvl, 1, j );

         for ( uint_t i = 1; i < n - 1 - j; ++i )
         {
            diagData[diagIdx] = stencils[dof][p1::stencil::C];
            ++diagIdx;
            ++dof;
         }
      }
   }

   void assemble_diagonalOperator_cell_precomputed_2d( std::shared_ptr< hyteg::Cell > cell, uint_t lvl, const real_t* diagData )
   {
      const auto& stencils = stencil_cell_3d_.at( cell->getID() )[lvl];

      uint_t dof = 0;
      for ( uint_t k = 1; k < n - 3; ++k )
      {
         for ( uint_t j = 1; j < n - 2 - k; ++j )
         {
            auto diagIdx = vertexdof::macrocell::index( lvl, 1, j, k );

            for ( uint_t i = 1; i < n - 1 - j - k; ++i )
            {
               diagData[diagIdx] = stencils[dof][p1::stencil::C];
               ++diagIdx;
               ++dof;
            }
         }
      }
   }

   void assemble_diagonalOperator_edge_surrogate_2d( std::shared_ptr< hyteg::Edge > edge, uint_t lvl, const real_t* diagData )
   {
      auto                       diagIdx   = vertexdof::macroedge::index( lvl, 1 );
      const PolyStencil< 2, 1 >& surrogate = surrogate_edge_2d_.at( edge->getID() )[lvl];
      const PolyDomain           X( lvl );
      for ( uint_t i = 1; i < n - 1; ++i )
      {
         const auto x = X[i];
         surrogate.eval( x );
         diagData[diagIdx] = surrogate.px()[p1::stencil::C];
         ++diagIdx;
      }
   }

   void assemble_diagonalOperator_face_surrogate_2d( std::shared_ptr< hyteg::Face > face, uint_t lvl, const real_t* diagData )
   {
      const PolyStencil< 2, 2 >& surrogate = surrogate_face_2d_.at( face->getID() )[lvl];
      const PolyDomain           X( lvl );

      for ( uint_t j = 1; j < n - 2; ++j )
      {
         auto diagIdx = vertexdof::macroface::index( lvl, 1, j );

         // restrict polynomial to 1D
         const auto y = X[j];
         surrogate.fix_y( y );

         for ( uint_t i = 1; i < n - 1 - j; ++i )
         {
            // evaluate polynomial
            const auto x = X[i];
            surrogate.eval( x );

            diagData[diagIdx] = surrogate.px()[p1::stencil::C];
            ++diagIdx;
         }
      }
   }

   void assemble_diagonalOperator_face_surrogate_3d( std::shared_ptr< hyteg::Face > face, uint_t lvl, const real_t* diagData )
   {
      const PolyStencil< 3, 2 >& surrogate = surrogate_face_3d_.at( face->getID() )[lvl];
      const PolyDomain           X( lvl );

      for ( uint_t j = 1; j < n - 2; ++j )
      {
         auto diagIdx = vertexdof::macroface::index( lvl, 1, j );

         // restrict polynomial to 1D
         const auto y = X[j];
         surrogate.fix_y( y );

         for ( uint_t i = 1; i < n - 1 - j; ++i )
         {
            // evaluate polynomial
            const auto x = X[i];
            surrogate.eval( x );

            diagData[diagIdx] = surrogate.px()[p1::stencil::C];
            ++diagIdx;
         }
      }
   }

   void assemble_diagonalOperator_cell_surrogate_2d( std::shared_ptr< hyteg::Cell > cell, uint_t lvl, const real_t* diagData )
   {
      const PolyStencil< 3, 3 >& surrogate = surrogate_cell_3d_.at( cell->getID() )[lvl];
      const PolyDomain           X( lvl );

      for ( uint_t k = 1; k < n - 3; ++k )
      {
         // restrict polynomial to 2D
         const auto z = X[k];
         surrogate.fix_z( z );

         for ( uint_t j = 1; j < n - 2 - k; ++j )
         {
            auto diagIdx = vertexdof::macrocell::index( lvl, 1, j, k );

            // restrict polynomial to 1D
            const auto y = X[j];
            surrogate.fix_y( y );

            for ( uint_t i = 1; i < n - 1 - j - k; ++i )
            {
               // evaluate polynomial
               const auto x = X[i];
               surrogate.eval( x );

               diagData[diagIdx] = surrogate.px()[p1::stencil::C];
               ++diagIdx;
            }
         }
      }
   }

   P1Form form_;

   bool is_initialized_;

   std::shared_ptr< P1Function< real_t > > diagonalValues_;
   std::shared_ptr< P1Function< real_t > > inverseDiagonalValues_;

   // least squares approximation for each level
   std::vector< std::shared_ptr< LSQ > > lsq_volume_;    // lsq for volume primitives (cells in 3d / faces in 2d)
   std::vector< std::shared_ptr< LSQ > > lsq_interface_; // lsq for interface primitives (faces in 3d / edges in 2d)
   uint_t                                downsampling_;

   // precomputed irregular stencils (vertices in 2D, vertices and edges in 3D
   VarStencilMap stencil_vtx_;
   VarStencilMap stencil_edge_3d_;
   // precomputed regular stencils for level 1-3
   StencilMap< 2 > stencil_edge_2d_;
   StencilMap< 2 > stencil_face_2d_;
   StencilMap< 3 > stencil_face_3d_;
   StencilMap< 3 > stencil_cell_3d_;
   // surrogates for level 4+
   PolyStencilMap< 2, 1 > surrogate_edge_2d_;
   PolyStencilMap< 2, 2 > surrogate_face_2d_;
   PolyStencilMap< 3, 2 > surrogate_face_3d_;
   PolyStencilMap< 3, 3 > surrogate_cell_3d_;
   // logical offsets on 3d faces, i.e. stencil directions as {-1,0,1}^3 tuples
   std::map< PrimitiveID, p1::stencil::StencilData< 3, indexing::Index > > offsets_face_3d_;
};

// todo test other forms

// typedef P1SurrogateOperator< forms::p1_diffusion_blending_q1 >  P1SurrogateLaplaceOperator;
// typedef P1SurrogateOperator< forms::p1_mass_blending_q4 >       P1SurrogateMassOperator;
// typedef P1SurrogateOperator< forms::p1_div_k_grad_blending_q3 > P1SurrogateDivkGradOperator;
// typedef P1SurrogateOperator< forms::p1_div_k_grad_affine_q3 >   P1SurrogateAffineDivkGradOperator;

} // namespace hyteg
