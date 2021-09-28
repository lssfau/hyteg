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
#pragma once

#include <array>

#include "hyteg/forms/P1LinearCombinationForm.hpp"
#include "hyteg/forms/P2LinearCombinationForm.hpp"
#include "hyteg/forms/P2RowSumForm.hpp"
#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/memory/LevelWiseMemory.hpp"
#include "hyteg/memory/StencilMemory.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFIndexing.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/sparseassembly/SparseMatrixProxy.hpp"
#include "hyteg/p1functionspace/P1Petsc.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"


namespace hyteg {

using walberla::real_t;

template < class P1Form, bool Diagonal = false, bool Lumped = false, bool InvertDiagonal = false >
class P1ConstantOperator : public Operator< P1Function< real_t >, P1Function< real_t > >,
                           public WeightedJacobiSmoothable< P1Function< real_t > >,
                           public GSSmoothable< P1Function< real_t > >,
                           public GSBackwardsSmoothable< P1Function< real_t > >,
                           public SORSmoothable< P1Function< real_t > >,
                           public SORBackwardsSmoothable< P1Function< real_t > >,
                           public OperatorWithInverseDiagonal< P1Function< real_t > >
{
 public:
   P1ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );
   P1ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel, const P1Form& form );

   ~P1ConstantOperator() override = default;

   void scale( real_t scalar );

   void apply( const P1Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const override final;

   void smooth_gs( const P1Function< real_t >& dst, const P1Function< real_t >& rhs, size_t level, DoFType flag ) const override;

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
                    const bool&                 backwards ) const;

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
                    const P1Function< real_t >& tmp,
                    real_t                      relax,
                    size_t                      level,
                    DoFType                     flag ) const override;

   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeDiagonalOperatorValues() { computeDiagonalOperatorValues( false ); }

   /// Trigger (re)computation of inverse diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   void computeInverseDiagonalOperatorValues() { computeDiagonalOperatorValues( true ); }

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

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P1Function< matIdx_t >&               src,
                  const P1Function< matIdx_t >&               dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      const auto storage = src.getStorage();

      for ( auto& it : this->getStorage()->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if ( testFlag( vertexBC, flag ) )
         {
            vertexdof::macrovertex::saveOperator(
                vertex, this->getVertexStencilID(), src.getVertexDataID(), dst.getVertexDataID(), mat, level );
         }
      }

      for ( auto& it : this->getStorage()->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            vertexdof::macroedge::saveOperator(
                level, edge, *storage, this->getEdgeStencilID(), src.getEdgeDataID(), dst.getEdgeDataID(), mat );
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
                  vertexdof::macroface::saveOperator3D(
                      level, face, *storage, this->getFaceStencil3DID(), src.getFaceDataID(), dst.getFaceDataID(), mat );
               }
               else
               {
                  vertexdof::macroface::saveOperator(
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
               vertexdof::macrocell::saveOperator(
                   level, cell, this->getCellStencilID(), src.getCellDataID(), dst.getCellDataID(), mat );
            }
         }
      }
   }

 private:
   void assembleStencils();

   void assembleStencils3D();

 private:
   /// Trigger (re)computation of diagonal matrix entries (central operator weights)
   /// Allocates the required memory if the function was not yet allocated.
   ///
   /// \param invert if true, assembles the function carrying the inverse of the diagonal
   void computeDiagonalOperatorValues( bool invert );

   std::shared_ptr< P1Function< real_t > > diagonalValues_;
   std::shared_ptr< P1Function< real_t > > inverseDiagonalValues_;

   void smooth_sor_macro_vertices( const P1Function< real_t >& dst,
                                   const P1Function< real_t >& rhs,
                                   real_t                      relax,
                                   size_t                      level,
                                   DoFType                     flag,
                                   const bool&                 backwards = false ) const;

   void smooth_sor_macro_edges( const P1Function< real_t >& dst,
                                const P1Function< real_t >& rhs,
                                real_t                      relax,
                                size_t                      level,
                                DoFType                     flag,
                                const bool&                 backwards = false ) const;

   void smooth_sor_macro_faces( const P1Function< real_t >& dst,
                                const P1Function< real_t >& rhs,
                                real_t                      relax,
                                size_t                      level,
                                DoFType                     flag,
                                const bool&                 backwards = false ) const;

   void smooth_sor_macro_cells( const P1Function< real_t >& dst,
                                const P1Function< real_t >& rhs,
                                real_t                      relax,
                                size_t                      level,
                                DoFType                     flag,
                                const bool&                 backwards = false ) const;

   PrimitiveDataID< StencilMemory< real_t >, Vertex >                             vertexStencilID_;
   PrimitiveDataID< StencilMemory< real_t >, Edge >                               edgeStencilID_;
   PrimitiveDataID< LevelWiseMemory< vertexdof::macroedge::StencilMap_T >, Edge > edgeStencil3DID_;
   PrimitiveDataID< StencilMemory< real_t >, Face >                               faceStencilID_;
   PrimitiveDataID< LevelWiseMemory< vertexdof::macroface::StencilMap_T >, Face > faceStencil3DID_;
   PrimitiveDataID< LevelWiseMemory< vertexdof::macrocell::StencilMap_T >, Cell > cellStencilID_;

   P1Form form_;
};

typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, fenics::NoAssemble > > P1ZeroOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise > >
    P1ConstantLaplaceOperator;
typedef P1ConstantOperator< P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, fenics::UndefinedAssembly >, true >
    P1DiagonalLaplaceOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_0_otherwise > > P1ConstantEpsilonOperator_11;
typedef P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_1_otherwise > > P1ConstantEpsilonOperator_12;
typedef P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_2_otherwise > > P1ConstantEpsilonOperator_21;
typedef P1ConstantOperator< P1FenicsForm< p1_stokes_epsilon_cell_integral_3_otherwise > > P1ConstantEpsilonOperator_22;

typedef P1ConstantOperator< P1FenicsForm< p1_div_cell_integral_0_otherwise, p1_tet_div_tet_cell_integral_0_otherwise > >
    P1DivxOperator;
typedef P1ConstantOperator< P1FenicsForm< p1_div_cell_integral_1_otherwise, p1_tet_div_tet_cell_integral_1_otherwise > >
                                                                                                           P1DivyOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_div_tet_cell_integral_2_otherwise > > P1DivzOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_divt_cell_integral_0_otherwise, p1_tet_divt_tet_cell_integral_0_otherwise > >
    P1DivTxOperator;
typedef P1ConstantOperator< P1FenicsForm< p1_divt_cell_integral_1_otherwise, p1_tet_divt_tet_cell_integral_1_otherwise > >
                                                                                                            P1DivTyOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_tet_divt_tet_cell_integral_2_otherwise > > P1DivTzOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise > >
    P1ConstantMassOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise >,
                            false,
                            true,
                            false >
    P1LumpedMassOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_mass_cell_integral_0_otherwise, p1_tet_mass_cell_integral_0_otherwise >,
                            false,
                            true,
                            true >
    P1LumpedInvMassOperator;

typedef P1ConstantOperator< P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise > >
    P1PSPGOperator;
typedef P1ConstantOperator< P1FenicsForm< p1_pspg_cell_integral_0_otherwise, p1_tet_pspg_tet_cell_integral_0_otherwise >,
                            true,
                            false,
                            true >
    P1PSPGInvDiagOperator;

typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_0_otherwise > >
    P2ToP1DivxVertexToVertexConstantOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p2_to_p1_tet_div_tet_cell_integral_1_otherwise > >
    P2ToP1DivyVertexToVertexConstantOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::UndefinedAssembly, p2_to_p1_tet_div_tet_cell_integral_2_otherwise > >
    P2ToP1DivzVertexToVertexConstantOperator;

typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_0_otherwise > >
    P1ToP1DivTxVertexToVertexConstantOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::NoAssemble, p1_to_p2_tet_divt_tet_cell_integral_1_otherwise > >
    P1ToP1DivTyVertexToVertexConstantOperator;
typedef P1ConstantOperator< P1FenicsForm< fenics::UndefinedAssembly, p1_to_p2_tet_divt_tet_cell_integral_2_otherwise > >
    P1ToP1DivTzVertexToVertexConstantOperator;

typedef P1ConstantOperator< P1LinearCombinationForm > P1ConstantLinearCombinationOperator;

} // namespace hyteg
