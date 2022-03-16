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
#include <hyteg/mixedoperators/polynomial/P2P1MacroFacePolynomial.hpp>
#include <hyteg/mixedoperators/polynomial/P2P1PolynomialDataHandling.hpp>
#include <hyteg/mixedoperators/variablestencil/P2P1VariableStencilCommon.hpp>
#include <hyteg/p1functionspace/VertexDoFFunction.hpp>
#include <hyteg/p2functionspace/P2Function.hpp>
#include <hyteg/p2functionspace/polynomial/StencilInterpolator.hpp>

#include "hyteg/forms/form_hyteg_generated/p2_to_p1/p2_to_p1_div_blending_q2.hpp"
#include "hyteg/types/pointnd.hpp"

namespace hyteg {

using walberla::real_t;

template < class P2ToP1Form >
class P2ToP1SurrogateOperator : public Operator< P2Function< real_t >, P1Function< real_t > >
{
 public:
   P2ToP1SurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                            size_t                                     minLevel,
                            size_t                                     maxLevel,
                            size_t                                     interpolationLevel )
   : Operator( storage, minLevel, maxLevel )
   , interpolationLevel_( interpolationLevel )
   {
      for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
      {
         PrimitiveDataID< P2toP1::FacePolynomialMemory, Face > id;
         auto dataHandling = std::make_shared< FaceP2toP1PolynomialMemoryDataHandling >();
         storage_->addFaceData( id, dataHandling, "P2toP1OperatorFacePolynomial" );
         polynomialIDs_[level] = id;
      }
   }

   P2ToP1SurrogateOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                            uint_t                                     minLevel,
                            uint_t                                     maxLevel,
                            uint_t                                     interpolationLevel,
                            uint_t                                     polyDegree )
   : P2ToP1SurrogateOperator( storage, minLevel, maxLevel, interpolationLevel )
   {
      interpolateStencils( polyDegree );
      useDegree( polyDegree );
   }

   ~P2ToP1SurrogateOperator() {}

   void interpolateStencils( uint_t polyDegree )
   {
      real_t H = 1.0 / ( walberla::real_c( levelinfo::num_microvertices_per_edge( interpolationLevel_ ) - 1 ) );

      // stencil entries
      std::array< real_t, P2::NumStencilentries2D::VtV > vertexToVertexStencil;
      std::array< real_t, P2::NumStencilentries2D::EtV > edgeToVertexStencil;

      // we only use polynomials for face stencils
      for ( auto& itF : storage_->getFaces() )
      {
         Face& face = *itF.second;
         form_.setGeometryMap( face.getGeometryMap() );

         Point3D x0( face.coords[0] ), x;

         Point3D D0 = face.coords[1] - face.coords[0];
         Point3D D2 = face.coords[2] - face.coords[0];

         for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
         {
            P2::StencilInterpolator< LSQPType::VERTEX, P2::NumStencilentries2D::VtV > VtVInterpolator( polyDegree,
                                                                                                       interpolationLevel_ );
            P2::StencilInterpolator< LSQPType::VERTEX, P2::NumStencilentries2D::EtV > EtVInterpolator( polyDegree,
                                                                                                       interpolationLevel_ );

            // directions (size of microfaces based on current level)
            real_t        h     = 1.0 / ( walberla::real_c( levelinfo::num_microvertices_per_edge( level ) - 1 ) );
            const Point3D dirS  = h * ( -D2 );
            const Point3D dirSE = h * ( D0 - D2 );
            const Point3D dirE  = h * ( D0 );
            const Point3D dirW  = -dirE;
            const Point3D dirNW = -dirSE;
            const Point3D dirN  = -dirS;

            // loop over all DOFs (number and position of microfaces based of interpolationLevel)
            for ( const auto& it : hyteg::edgedof::macroface::Iterator( interpolationLevel_, 0 ) )
            {
               // position of vertex on macroface
               real_t r = walberla::real_c( it.row() ) * H;
               real_t c = walberla::real_c( it.col() ) * H;
               x        = x0 + ( r * D2 + c * D0 );
               // corresponding point on reference element
               Point2D xi( {c, r} );

               P2toP1::variablestencil::macroface::assembleStencil(
                   form_, x, dirS, dirSE, dirE, dirN, dirNW, dirW, vertexToVertexStencil, edgeToVertexStencil );

               // add interpolation points
               // vertex DoF
               if ( !vertexdof::macroface::isVertexOnBoundary( interpolationLevel_, it ) )
               {
                  for ( uint_t i = 0; i < P2::NumStencilentries2D::VtV; ++i )
                  {
                     VtVInterpolator[i].addInterpolationPoint( xi, vertexToVertexStencil[i] );
                  }

                  for ( uint_t i = 0; i < P2::NumStencilentries2D::EtV; ++i )
                  {
                     EtVInterpolator[i].addInterpolationPoint( xi, edgeToVertexStencil[i] );
                  }
               }
            }

            auto id = face.getData( polynomialIDs_[level] );

            auto& poly = id->addDegree( polyDegree );

            VtVInterpolator.interpolate( poly.VtV );
            EtVInterpolator.interpolate( poly.EtV );
         }
      }
   }

   void useDegree( uint_t degree ) { polyDegree_ = degree; }

   void apply( const P2Function< real_t >& src,
               const P1Function< real_t >& dst,
               size_t                      level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const
   {
      checkForMissingPolynomial( level );

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
            P2toP1::variablestencil::macroface::applyPolynomial( polyDegree_,
                                                                 polynomialIDs_.at( level ),
                                                                 level,
                                                                 face,
                                                                 srcVertexDoF.getFaceDataID(),
                                                                 srcEdgeDoF.getFaceDataID(),
                                                                 dst.getFaceDataID(),
                                                                 updateType );
         }
      }

      if ( storage_->hasGlobalCells() )
      {
         WALBERLA_ABORT( "P2ToP1SurrogateOperator not implemented for 3D" )
      }
   }

   void smooth_gs( P2Function< real_t >& dst, P1Function< real_t >& rhs, size_t level, DoFType flag )
   {
      WALBERLA_ABORT( "not implemented" );
   }

 private:
   inline void checkForMissingPolynomial( uint_t level ) const
   {
      WALBERLA_ASSERT( polynomialIDs_.count( level ) > 0, "Polynomial for level " << level << " has not been interpolated" );
   }

   uint_t                                                                    polyDegree_;
   uint_t                                                                    interpolationLevel_;
   P2ToP1Form                                                                form_;
   std::map< uint_t, PrimitiveDataID< P2toP1::FacePolynomialMemory, Face > > polynomialIDs_;
};

typedef P2ToP1SurrogateOperator< forms::p2_to_p1_div_0_blending_q2 > P2ToP1SurrogateDivxOperator;
typedef P2ToP1SurrogateOperator< forms::p2_to_p1_div_1_blending_q2 > P2ToP1SurrogateDivyOperator;
typedef P2ToP1SurrogateOperator< forms::p2_to_p1_div_2_blending_q2 > P2ToP1SurrogateDivzOperator;

} // namespace hyteg
