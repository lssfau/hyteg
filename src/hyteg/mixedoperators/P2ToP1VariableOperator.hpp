/*
 * Copyright (c) 2017-2019 Benjamin Mann.
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
#include <hyteg/communication/Syncing.hpp>
#include <hyteg/mixedoperators/variablestencil/P2P1VariableStencilCommon.hpp>
#include <hyteg/p1functionspace/VertexDoFFunction.hpp>
#include <hyteg/p2functionspace/P2Function.hpp>

#include "hyteg/forms/form_hyteg_generated/p2_to_p1/p2_to_p1_div_blending_q2.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

using walberla::real_t;

template < class P2ToP1Form >
class P2ToP1VariableOperator : public Operator< P2Function< real_t >, P1Function< real_t > >
{
 public:
   P2ToP1VariableOperator( const std::shared_ptr< PrimitiveStorage >& storage, size_t minLevel, size_t maxLevel )
   : Operator( storage, minLevel, maxLevel )
   {}

   void apply( const P2Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      communication::syncP2FunctionBetweenPrimitives( src, level );

      const vertexdof::VertexDoFFunction< real_t >& srcVertexDoF = src.getVertexDoFFunction();
      const EdgeDoFFunction< real_t >&              srcEdgeDoF   = src.getEdgeDoFFunction();

      for ( auto& it : storage_->getVertices() )
      {
         hyteg::Vertex& vertex = *it.second;

         const DoFType vtxFlag = dst.getBoundaryCondition().getBoundaryType( vertex.getMeshBoundaryFlag() );

         if ( testFlag( vtxFlag, flag ) )
         {
            P2toP1::variablestencil::macrovertex::applyVariableStencil< P2ToP1Form >( level,
                                                                                      vertex,
                                                                                      storage_,
                                                                                      srcVertexDoF.getVertexDataID(),
                                                                                      srcEdgeDoF.getVertexDataID(),
                                                                                      dst.getVertexDataID(),
                                                                                      updateType );
         }
      }

      for ( auto& it : storage_->getEdges() )
      {
         hyteg::Edge& edge = *it.second;

         const DoFType edgeFlag = dst.getBoundaryCondition().getBoundaryType( edge.getMeshBoundaryFlag() );

         if ( testFlag( edgeFlag, flag ) )
         {
            P2toP1::variablestencil::macroedge::applyVariableStencil< P2ToP1Form >( level,
                                                                                    edge,
                                                                                    storage_,
                                                                                    srcVertexDoF.getEdgeDataID(),
                                                                                    srcEdgeDoF.getEdgeDataID(),
                                                                                    dst.getEdgeDataID(),
                                                                                    updateType );
         }
      }

      for ( auto& it : storage_->getFaces() )
      {
         hyteg::Face& face = *it.second;

         const DoFType faceFlag = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );

         if ( testFlag( faceFlag, flag ) )
         {
            P2toP1::variablestencil::macroface::applyVariableStencil< P2ToP1Form >(
                level, face, srcVertexDoF.getFaceDataID(), srcEdgeDoF.getFaceDataID(), dst.getFaceDataID(), updateType );
         }
      }

      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "P2ToP1VariableOperator not implemented for 3D" )
      }
   }

   void smooth_gs( P2Function< real_t >& dst, P1Function< real_t >& rhs, size_t level, DoFType flag )
   {
      WALBERLA_ABORT( "not implemented" );
   }
};

typedef P2ToP1VariableOperator< forms::p2_to_p1_div_0_blending_q2 > P2ToP1BlendingDivxOperator;
typedef P2ToP1VariableOperator< forms::p2_to_p1_div_1_blending_q2 > P2ToP1BlendingDivyOperator;
typedef P2ToP1VariableOperator< forms::p2_to_p1_div_2_blending_q2 > P2ToP1BlendingDivzOperator;

} // namespace hyteg
