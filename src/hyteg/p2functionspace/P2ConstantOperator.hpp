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
#pragma once

#include "hyteg/edgedofspace/EdgeDoFOperator.hpp"
#include "hyteg/forms/P2LinearCombinationForm.hpp"
#include "hyteg/forms/P2RowSumForm.hpp"
#include "hyteg/forms/form_fenics_base/P2FenicsForm.hpp"
#include "hyteg/forms/P1WrapperForm.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFOperator.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/VertexDoFToEdgeDoFOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/solvers/Smoothables.hpp"

namespace hyteg {

using walberla::real_t;

template < class P2Form >
class P2ConstantOperator : public Operator< P2Function< real_t >, P2Function< real_t > >,
                           public WeightedJacobiSmoothable< P2Function< real_t > >,
                           public GSSmoothable< P2Function< real_t > >,
                           public GSBackwardsSmoothable< P2Function< real_t > >,
                           public SORSmoothable< P2Function< real_t > >,
                           public SORBackwardsSmoothable< P2Function< real_t > >
{
 public:
   P2ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel );
   P2ConstantOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel, const P2Form& form );

   const P1ConstantOperator< P1WrapperForm<P2Form> >& getVertexToVertexOpr() const { return vertexToVertex; }

   const EdgeDoFToVertexDoFOperator< P2Form >& getEdgeToVertexOpr() const { return edgeToVertex; }

   const VertexDoFToEdgeDoFOperator< P2Form >& getVertexToEdgeOpr() const { return vertexToEdge; }

   const EdgeDoFOperator< P2Form >& getEdgeToEdgeOpr() const { return edgeToEdge; }

   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const override final;

   void smooth_gs( const P2Function< real_t >& dst, const P2Function< real_t >& rhs, size_t level, DoFType flag ) const override;

   void smooth_gs_backwards( const P2Function< real_t >& dst,
                             const P2Function< real_t >& rhs,
                             size_t                      level,
                             DoFType                     flag ) const override
   {
      if ( !storage_->hasGlobalCells() )
      {
         throw std::runtime_error( "P2ConstantOperator: Backward GS currently only implemented for 3D." );
      }
      smooth_sor_backwards( dst, rhs, 1.0, level, flag );
   }

   void smooth_sor( const P2Function< real_t >& dst,
                    const P2Function< real_t >& rhs,
                    real_t                      relax,
                    size_t                      level,
                    DoFType                     flag,
                    const bool&                 backwards ) const;

   void smooth_sor( const P2Function< real_t >& dst,
                    const P2Function< real_t >& rhs,
                    real_t                      relax,
                    size_t                      level,
                    DoFType                     flag ) const override
   {
      smooth_sor( dst, rhs, relax, level, flag, false );
   }

   void smooth_sor_backwards( const P2Function< real_t >& dst,
                              const P2Function< real_t >& rhs,
                              real_t                      relax,
                              size_t                      level,
                              DoFType                     flag ) const override
   {
      if ( !storage_->hasGlobalCells() )
      {
         throw std::runtime_error( "P2ConstantOperator: Backward SOR currently only implemented for 3D." );
      }
      smooth_sor( dst, rhs, relax, level, flag, true );
   }

   void smooth_jac( const P2Function< real_t >& dst,
                    const P2Function< real_t >& rhs,
                    const P2Function< real_t >& src,
                    real_t                      relax,
                    size_t                      level,
                    DoFType                     flag ) const override;

   void toMatrix( const std::shared_ptr< SparseMatrixProxy >& mat,
                  const P2Function< idx_t >&                  src,
                  const P2Function< idx_t >&                  dst,
                  size_t                                      level,
                  DoFType                                     flag ) const override
   {
      this->getVertexToVertexOpr().toMatrix( mat, src.getVertexDoFFunction(), dst.getVertexDoFFunction(), level, flag );
      this->getEdgeToVertexOpr().toMatrix( mat, src.getEdgeDoFFunction(), dst.getVertexDoFFunction(), level, flag );
      this->getVertexToEdgeOpr().toMatrix( mat, src.getVertexDoFFunction(), dst.getEdgeDoFFunction(), level, flag );
      this->getEdgeToEdgeOpr().toMatrix( mat, src.getEdgeDoFFunction(), dst.getEdgeDoFFunction(), level, flag );
   }

 private:
   void smooth_sor_macro_vertices( const P2Function< real_t >& dst,
                                   const P2Function< real_t >& rhs,
                                   const real_t&               relax,
                                   size_t                      level,
                                   DoFType                     flag,
                                   const bool&                 backwards = false ) const;

   void smooth_sor_macro_edges( const P2Function< real_t >& dst,
                                const P2Function< real_t >& rhs,
                                const real_t&               relax,
                                size_t                      level,
                                DoFType                     flag,
                                const bool&                 backwards = false ) const;

   void smooth_sor_macro_faces( const P2Function< real_t >& dst,
                                const P2Function< real_t >& rhs,
                                const real_t&               relax,
                                size_t                      level,
                                DoFType                     flag,
                                const bool&                 backwards = false ) const;

   void smooth_sor_macro_cells( const P2Function< real_t >& dst,
                                const P2Function< real_t >& rhs,
                                const real_t&               relax,
                                size_t                      level,
                                DoFType                     flag,
                                const bool&                 backwards = false ) const;

   P1ConstantOperator< P1WrapperForm<P2Form> >         vertexToVertex;
   EdgeDoFToVertexDoFOperator< P2Form > edgeToVertex;
   VertexDoFToEdgeDoFOperator< P2Form > vertexToEdge;
   EdgeDoFOperator< P2Form >            edgeToEdge;

   P2Form form_;
};

typedef P2ConstantOperator< P2FenicsForm< p2_diffusion_cell_integral_0_otherwise, p2_tet_diffusion_cell_integral_0_otherwise > >
    P2ConstantLaplaceOperator;
typedef P2ConstantOperator< P2FenicsForm< p2_mass_cell_integral_0_otherwise, p2_tet_mass_cell_integral_0_otherwise > >
    P2ConstantMassOperator;

typedef P2ConstantOperator< P2FenicsForm< p2_divt_cell_integral_0_otherwise, p2_tet_divt_tet_cell_integral_0_otherwise > >
    P2ConstantDivTxOperator;
typedef P2ConstantOperator< P2FenicsForm< p2_divt_cell_integral_1_otherwise, p2_tet_divt_tet_cell_integral_1_otherwise > >
    P2ConstantDivTyOperator;
typedef P2ConstantOperator< P2FenicsForm< fenics::NoAssemble, p2_tet_divt_tet_cell_integral_2_otherwise > >
    P2ConstantDivTzOperator;
typedef P2ConstantOperator< P2FenicsForm< p2_div_cell_integral_0_otherwise, p2_tet_div_tet_cell_integral_0_otherwise > >
    P2ConstantDivxOperator;
typedef P2ConstantOperator< P2FenicsForm< p2_div_cell_integral_1_otherwise, p2_tet_div_tet_cell_integral_1_otherwise > >
                                                                                                           P2ConstantDivyOperator;
typedef P2ConstantOperator< P2FenicsForm< fenics::NoAssemble, p2_tet_div_tet_cell_integral_2_otherwise > > P2ConstantDivzOperator;

typedef P2ConstantOperator< P2FenicsForm< p2_pspg_cell_integral_0_otherwise, p2_tet_pspg_tet_cell_integral_0_otherwise > >
    P2ConstantPSPGOperator;

typedef P2ConstantOperator< P2LinearCombinationForm > P2ConstantLinearCombinationOperator;
typedef P2ConstantOperator< P2RowSumForm >            P2ConstantRowSumOperator;

} // namespace hyteg
