/*
 * Copyright (c) 2017-2023 Daniel Drzisga, Dominik Thoennes, Marcus Mohr.
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_mass_blending_q4.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormDivKGrad.hpp"
#include "hyteg/forms/form_hyteg_manual/P2FormLaplace.hpp"
#include "hyteg/operators/Operator.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/variablestencil/P2VariableStencilCommon.hpp"
#include "hyteg/types/PointND.hpp"
#include "hyteg/solvers/Smoothables.hpp"

namespace hyteg {

template < class P2Form >
class P2VariableOperator : public Operator< P2Function< real_t >, P2Function< real_t > >,
                           public GSSmoothable< P2Function< real_t > >,
                           public ConstantJacobiSmoothable< P2Function< real_t > >
{
 public:
   P2VariableOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   {}

   ~P2VariableOperator() override = default;

   void apply( const P2Function< real_t >& src,
               const P2Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

      communication::syncFunctionBetweenPrimitives( src, level );

      const vertexdof::VertexDoFFunction< real_t >& srcVertexDoF = src.getVertexDoFFunction();
      const EdgeDoFFunction< real_t >&              srcEdgeDoF   = src.getEdgeDoFFunction();
      const vertexdof::VertexDoFFunction< real_t >& dstVertexDoF = dst.getVertexDoFFunction();
      const EdgeDoFFunction< real_t >&              dstEdgeDoF   = dst.getEdgeDoFFunction();

      for ( auto& it : storage_->getVertices() )
      {
         hyteg::Vertex& vertex = *it.second;

         const DoFType vtxFlag = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );

         if ( testFlag( vtxFlag, flag ) )
         {
            P2::variablestencil::macrovertex::applyVariableStencil< P2Form >( level,
                                                                              vertex,
                                                                              storage_,
                                                                              srcVertexDoF.getVertexDataID(),
                                                                              srcEdgeDoF.getVertexDataID(),
                                                                              dstVertexDoF.getVertexDataID(),
                                                                              updateType );
         }
      }

      for ( auto& it : storage_->getEdges() )
      {
         hyteg::Edge& edge = *it.second;

         const DoFType edgeFlag = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );

         if ( testFlag( edgeFlag, flag ) )
         {
            P2::variablestencil::macroedge::applyVariableStencil< P2Form >( level,
                                                                            edge,
                                                                            storage_,
                                                                            srcVertexDoF.getEdgeDataID(),
                                                                            srcEdgeDoF.getEdgeDataID(),
                                                                            dstVertexDoF.getEdgeDataID(),
                                                                            dstEdgeDoF.getEdgeDataID(),
                                                                            updateType );
         }
      }

      for ( auto& it : storage_->getFaces() )
      {
         hyteg::Face& face = *it.second;

         const DoFType faceFlag = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );

         if ( testFlag( faceFlag, flag ) )
         {
            P2::variablestencil::macroface::applyVariableStencil< P2Form >( level,
                                                                            face,
                                                                            srcVertexDoF.getFaceDataID(),
                                                                            srcEdgeDoF.getFaceDataID(),
                                                                            dstVertexDoF.getFaceDataID(),
                                                                            dstEdgeDoF.getFaceDataID(),
                                                                            updateType );
         }
      }
   }

   void smooth_gs( const P2Function< real_t >& dst, const P2Function< real_t >& rhs, size_t level, DoFType flag ) const override
   {
      communication::syncFunctionBetweenPrimitives( dst, level );

      const vertexdof::VertexDoFFunction< real_t >& dstVertexDoF = dst.getVertexDoFFunction();
      const EdgeDoFFunction< real_t >&              dstEdgeDoF   = dst.getEdgeDoFFunction();
      const vertexdof::VertexDoFFunction< real_t >& rhsVertexDoF = rhs.getVertexDoFFunction();
      const EdgeDoFFunction< real_t >&              rhsEdgeDoF   = rhs.getEdgeDoFFunction();

      for ( auto& it : storage_->getVertices() )
      {
         hyteg::Vertex& vertex = *it.second;

         const DoFType vertexFlag = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );

         if ( testFlag( vertexFlag, flag ) )
         {
            P2::variablestencil::macrovertex::smoothGSVariableStencil< P2Form >( level,
                                                                                 vertex,
                                                                                 storage_,
                                                                                 dstVertexDoF.getVertexDataID(),
                                                                                 dstEdgeDoF.getVertexDataID(),
                                                                                 rhsVertexDoF.getVertexDataID() );
         }
      }

      communication::syncFunctionBetweenPrimitives( dst, level );

      for ( auto& it : storage_->getEdges() )
      {
         hyteg::Edge& edge = *it.second;

         const DoFType edgeFlag = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );

         if ( testFlag( edgeFlag, flag ) )
         {
            P2::variablestencil::macroedge::smoothGSVariableStencil< P2Form >( level,
                                                                               edge,
                                                                               storage_,
                                                                               dstVertexDoF.getEdgeDataID(),
                                                                               dstEdgeDoF.getEdgeDataID(),
                                                                               rhsVertexDoF.getEdgeDataID(),
                                                                               rhsEdgeDoF.getEdgeDataID() );
         }
      }

      communication::syncFunctionBetweenPrimitives( dst, level );

      for ( auto& it : storage_->getFaces() )
      {
         hyteg::Face& face = *it.second;

         const DoFType faceFlag = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );

         if ( testFlag( faceFlag, flag ) )
         {
            P2::variablestencil::macroface::smoothGSVariableStencil< P2Form >( level,
                                                                               face,
                                                                               dstVertexDoF.getFaceDataID(),
                                                                               dstEdgeDoF.getFaceDataID(),
                                                                               rhsVertexDoF.getFaceDataID(),
                                                                               rhsEdgeDoF.getFaceDataID() );
         }
      }

      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "P2VariableOperator not implemented for 3D" )
      }
   }

   void smooth_jac( const P2Function< real_t >& dst,
                    const P2Function< real_t >& rhs,
                    const P2Function< real_t >& tmp,
                    size_t                      level,
                    DoFType                     flag ) const override
   {
      WALBERLA_ABORT( "To be implemented" );
   }
};
typedef P2VariableOperator< P2Form_laplace >             P2BlendingLaplaceOperator;
typedef P2VariableOperator< P2Form_divKgrad >            P2divKgradOperator;
typedef P2VariableOperator< forms::p2_mass_blending_q4 > P2BlendingMassOperator;

} // namespace hyteg
