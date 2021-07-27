/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr.
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
#include <hyteg/operators/Operator.hpp>
#include <hyteg/p1functionspace/VertexDoFMacroEdge.hpp>
#include <hyteg/p1functionspace/VertexDoFMacroFace.hpp>
#include <hyteg/p1functionspace/VertexDoFMacroVertex.hpp>
#include <hyteg/p1functionspace/variablestencil/VertexDoFVariableStencil.hpp>

#include "hyteg/forms/form_hyteg_generated/P1FormDiv.hpp"
#include "hyteg/forms/form_hyteg_generated/P1FormDivT.hpp"
#include "hyteg/forms/form_hyteg_generated/deprecated/P1FormPSPG.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_epsilon_all_forms.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_mass_blending_q4.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/solvers/Smoothables.hpp"
#include "hyteg/types/pointnd.hpp"

#include "P1DataHandling.hpp"

namespace hyteg {

template < class P1Form >
class P1VariableOperator : public Operator< P1Function< real_t >, P1Function< real_t > >,
                           public GSSmoothable< P1Function< real_t > >,
                           public ConstantJacobiSmoothable< P1Function< real_t > >
{
 public:
   P1VariableOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   {}

   ~P1VariableOperator() override = default;

   void apply( const P1Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const override final
   {
      WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

      src.communicate< Vertex, Edge >( level );
      src.communicate< Edge, Face >( level );
      src.communicate< Face, Edge >( level );
      src.communicate< Edge, Vertex >( level );

      for ( auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if ( testFlag( vertexBC, flag ) )
         {
            vertexdof::variablestencil::macrovertex::applyVariableStencil< real_t, P1Form >(
                level, vertex, storage_, src.getVertexDataID(), dst.getVertexDataID(), updateType );
         }
      }

      for ( auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            vertexdof::variablestencil::macroedge::applyVariableStencil< real_t, P1Form >(
                level, edge, storage_, src.getEdgeDataID(), dst.getEdgeDataID(), updateType );
         }
      }

      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            vertexdof::variablestencil::macroface::applyVariableStencil< real_t, P1Form >(
                level, face, src.getFaceDataID(), dst.getFaceDataID(), updateType );
         }
      }
   }

   void smooth_gs( const P1Function< real_t >& dst, const P1Function< real_t >& rhs, size_t level, DoFType flag ) const override
   {
      dst.communicate< Vertex, Edge >( level );
      dst.communicate< Edge, Face >( level );

      dst.communicate< Face, Edge >( level );
      dst.communicate< Edge, Vertex >( level );

      for ( auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if ( testFlag( vertexBC, flag ) )
         {
            vertexdof::variablestencil::macrovertex::smoothGSVariableStencil< real_t, P1Form >(
                level, vertex, storage_, dst.getVertexDataID(), rhs.getVertexDataID() );
         }
      }

      dst.communicate< Vertex, Edge >( level );

      for ( auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            vertexdof::variablestencil::macroedge::smoothGSVariableStencil< real_t, P1Form >(
                level, edge, storage_, dst.getEdgeDataID(), rhs.getEdgeDataID() );
         }
      }

      dst.communicate< Edge, Face >( level );

      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            vertexdof::variablestencil::macroface::smoothGSVariableStencil< real_t, P1Form >(
                level, face, dst.getFaceDataID(), rhs.getFaceDataID() );
         }
      }
      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "P1VariableOperator not implemented for 3D" )
      }
   }

   void smooth_jac( const P1Function< real_t >& dst,
                    const P1Function< real_t >& rhs,
                    const P1Function< real_t >& tmp,
                    size_t                      level,
                    DoFType                     flag ) const override
   {
      // start pulling vertex halos
      tmp.startCommunication< Edge, Vertex >( level );

      // start pulling edge halos
      tmp.startCommunication< Face, Edge >( level );

      // end pulling vertex halos
      tmp.endCommunication< Edge, Vertex >( level );

      for ( auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if ( testFlag( vertexBC, flag ) )
         {
            WALBERLA_ABORT( "To be implemented" )
            //        P1Vertex::smooth_jac(vertex, vertexLocalMatrixID_, dst.getVertexDataID(), rhs.getVertexDataID(), tmp.getVertexDataID(), level);
         }
      }

      dst.startCommunication< Vertex, Edge >( level );

      // end pulling edge halos
      tmp.endCommunication< Face, Edge >( level );

      for ( auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            WALBERLA_ABORT( "To be implemented" )
            //        P1Edge::smooth_jac(level, edge, edgeLocalMatrixID_, dst.getEdgeDataID(), rhs.getEdgeDataID(), tmp.getEdgeDataID());
         }
      }

      dst.endCommunication< Vertex, Edge >( level );

      dst.startCommunication< Edge, Face >( level );

      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            WALBERLA_ABORT( "To be implemented" )
            //        P1Face::smooth_jac(level, face, faceLocalMatrixID_, dst.getFaceDataID(), rhs.getFaceDataID(), tmp.getFaceDataID());
         }
      }

      dst.endCommunication< Edge, Face >( level );
   }

#ifdef HYTEG_BUILD_WITH_PETSC
   void createMatrix_impl( P1Function< real_t >& src, P1Function< real_t >& dst, Mat& mat, size_t level, DoFType flag )
   {
      for ( auto& it : storage_->getVertices() )
      {
         Vertex& vertex = *it.second;

         const DoFType vertexBC = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );
         if ( testFlag( vertexBC, flag ) )
         {
            WALBERLA_ABORT( "To be implemented" )
            //        P1Vertex::saveOperator(vertex, vertexLocalMatrixID_, src.getVertexDataID(), dst.getVertexDataID(), mat, level);
         }
      }

      for ( auto& it : storage_->getEdges() )
      {
         Edge& edge = *it.second;

         const DoFType edgeBC = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );
         if ( testFlag( edgeBC, flag ) )
         {
            WALBERLA_ABORT( "To be implemented" )
            //        P1Edge::saveOperator(level, edge, edgeLocalMatrixID_, src.getEdgeDataID(), dst.getEdgeDataID(), mat);
         }
      }

      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            WALBERLA_ABORT( "To be implemented" )
            //        P1Face::saveOperator(level, face, faceLocalMatrixID_, src.getFaceDataID(), dst.getFaceDataID(), mat);
         }
      }
   }
#endif
};

typedef P1VariableOperator< forms::p1_diffusion_blending_q3 > P1BlendingLaplaceOperator;
typedef P1VariableOperator< forms::p1_mass_blending_q4 >      P1BlendingMassOperator;

typedef P1VariableOperator< forms::p1_epsiloncc_0_0_blending_q2 > P1BlendingEpsilonOperator_11;
typedef P1VariableOperator< forms::p1_epsiloncc_0_1_blending_q2 > P1BlendingEpsilonOperator_12;
typedef P1VariableOperator< forms::p1_epsiloncc_1_0_blending_q2 > P1BlendingEpsilonOperator_21;
typedef P1VariableOperator< forms::p1_epsiloncc_1_1_blending_q2 > P1BlendingEpsilonOperator_22;

typedef P1VariableOperator< P1Form_divT_1 > P1BlendingDivTOperator_1;
typedef P1VariableOperator< P1Form_divT_2 > P1BlendingDivTOperator_2;

typedef P1VariableOperator< P1Form_div_1 > P1BlendingDivOperator_1;
typedef P1VariableOperator< P1Form_div_2 > P1BlendingDivOperator_2;

typedef P1VariableOperator< P1Form_pspg > P1BlendingPSPGOperator;

} // namespace hyteg
