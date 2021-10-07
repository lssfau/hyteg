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
#include <hyteg/operators/Operator.hpp>

#include "hyteg/types/pointnd.hpp"

#include "P1DataHandling.hpp"

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

#include "variablestencil/VertexDoFVariableStencil.hpp"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

#include "hyteg/forms/form_hyteg_generated/p1/p1_div_0_blending_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_1_blending_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_divt_0_blending_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_divt_1_blending_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_epsilon_all_forms.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_diffusion_blending_q1.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_mass_blending_q4.hpp"
#include "hyteg/forms/form_hyteg_generated/deprecated/P1FormPSPG.hpp"
#include "hyteg/p1functionspace/VertexDoFMemory.hpp"
#include "hyteg/polynomial/LSQPInterpolator.hpp"

#include "VertexDoFMacroEdge.hpp"
#include "VertexDoFMacroFace.hpp"
#include "VertexDoFMacroVertex.hpp"
#include "polynomial/VertexDoFMacroFacePolynomial.hpp"

namespace hyteg {

template < class P1Form, OperatorType OprType >
class P1PolynomialBlendingOperator : public Operator< P1Function< real_t >, P1Function< real_t > >,
                                     public GSSmoothable< P1Function< real_t > >,
                                     public ConstantJacobiSmoothable< P1Function< real_t > >
{
 public:
   typedef LSQPInterpolator< MonomialBasis2D, LSQPType::EDGE >   EdgeInterpolator;
   typedef LSQPInterpolator< MonomialBasis2D, LSQPType::VERTEX > VertexInterpolator;

   P1PolynomialBlendingOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                 uint_t                                     minLevel,
                                 uint_t                                     maxLevel,
                                 uint_t                                     interpolationLevel )
   : Operator( storage, minLevel, maxLevel )
   , interpolationLevel_( interpolationLevel )
   {
      for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
      {
         PrimitiveDataID< FaceP1PolynomialMemory, Face > facePolynomialID;
         auto faceP1PolynomialMemoryDataHandling = std::make_shared< FaceP1PolynomialMemoryDataHandling >( polyDegree_ );
         storage_->addFaceData( facePolynomialID, faceP1PolynomialMemoryDataHandling, "P1OperatorFacePolynomial" );
         facePolynomialIDs_[level] = facePolynomialID;
      }
   }

   P1PolynomialBlendingOperator( const std::shared_ptr< PrimitiveStorage >& storage,
                                 uint_t                                     minLevel,
                                 uint_t                                     maxLevel,
                                 uint_t                                     interpolationLevel,
                                 uint_t                                     polyDegree )
   : P1PolynomialBlendingOperator( storage, minLevel, maxLevel, interpolationLevel )
   {
      interpolateStencils( polyDegree );
      useDegree( polyDegree );
   }

   ~P1PolynomialBlendingOperator() {}

   void interpolateStencils( uint_t polyDegree )
   {
      if ( OprType == OperatorType::ODD )
      {
         interpolateStencilsAsymmetric( polyDegree );
      }
      else
      {
         interpolateStencilsSymmetric( polyDegree );
      }
   }

   void interpolateStencilsSymmetric( uint_t polyDegree )
   {
      typedef stencilDirection SD;
      using namespace P1Elements;

      std::vector< real_t > faceStencil( 7 );

      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;
         form.setGeometryMap( face.getGeometryMap() );

         for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
         {
            auto facePolynomials = face.getData( facePolynomialIDs_[level] );
            facePolynomials->addDegree( polyDegree );

            uint_t rowsize       = levelinfo::num_microvertices_per_edge( interpolationLevel_ );
            uint_t rowsizeFine   = levelinfo::num_microvertices_per_edge( level );
            uint_t inner_rowsize = rowsize;

            VertexInterpolator centerInterpolator( polyDegree, interpolationLevel_ );
            EdgeInterpolator   horiInterpolator( polyDegree, interpolationLevel_ );
            EdgeInterpolator   vertInterpolator( polyDegree, interpolationLevel_ );
            EdgeInterpolator   diagInterpolator( polyDegree, interpolationLevel_ );

            Point3D x, x0;
            x0 = face.coords[0];

            real_t ref_H = 1.0 / ( walberla::real_c( rowsize - 1 ) );
            real_t ref_h = 1.0 / ( walberla::real_c( rowsizeFine - 1 ) );

            Point3D d0 = ref_H * ( face.coords[1] - face.coords[0] );
            Point3D d2 = ref_H * ( face.coords[2] - face.coords[0] );

            // fine directions
            Point3D d0f = ref_h * ( face.coords[1] - face.coords[0] );
            Point3D d2f = ref_h * ( face.coords[2] - face.coords[0] );

            Point2D ref_x;

            Point3D dirS  = -1.0 * d2f;
            Point3D dirSE = d0f - 1.0 * d2f;
            Point3D dirE  = d0f;
            Point3D dirW  = -1.0 * d0f;
            Point3D dirNW = -1.0 * d0f + d2f;
            Point3D dirN  = d2f;

            for ( uint_t j = 1; j < rowsize - 2; ++j )
            {
               ref_x[1] = j * ref_H;

               x = x0;
               x += walberla::real_c( j ) * d2 + d0;

               uint_t i;
               for ( i = 1; i < inner_rowsize - 2; ++i )
               {
                  ref_x[0] = i * ref_H;

                  std::fill( faceStencil.begin(), faceStencil.end(), 0.0 );

                  vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                      form, { x, x + dirW, x + dirS }, P1Elements::P1Elements2D::elementSW, faceStencil.data() );
                  vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                      form, { x, x + dirS, x + dirSE }, P1Elements::P1Elements2D::elementS, faceStencil.data() );
                  vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                      form, { x, x + dirSE, x + dirE }, P1Elements::P1Elements2D::elementSE, faceStencil.data() );
                  vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                      form, { x, x + dirE, x + dirN }, P1Elements::P1Elements2D::elementNE, faceStencil.data() );
                  vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                      form, { x, x + dirN, x + dirNW }, P1Elements::P1Elements2D::elementN, faceStencil.data() );
                  vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                      form, { x, x + dirNW, x + dirW }, P1Elements::P1Elements2D::elementNW, faceStencil.data() );

                  //if (i == 1 && j == 1) {
                  //   PointND<real_t, 7> test(faceStencil.data());
                  //   WALBERLA_LOG_INFO("stencil = " << test);
                  //}

                  centerInterpolator.addInterpolationPoint( ref_x,
                                                            faceStencil[vertexdof::stencilIndexFromVertex( SD::VERTEX_C )] );

                  horiInterpolator.addInterpolationPoint( ref_x + Point2D{ { 0.5 * ref_h, 0.0 } },
                                                          faceStencil[vertexdof::stencilIndexFromVertex( SD::VERTEX_E )] );

                  vertInterpolator.addInterpolationPoint( ref_x + Point2D{ { 0.0, -0.5 * ref_h } },
                                                          faceStencil[vertexdof::stencilIndexFromVertex( SD::VERTEX_S )] );

                  diagInterpolator.addInterpolationPoint( ref_x + Point2D{ { 0.5 * ref_h, -0.5 * ref_h } },
                                                          faceStencil[vertexdof::stencilIndexFromVertex( SD::VERTEX_SE )] );

                  if ( i == 1 )
                  {
                     diagInterpolator.addInterpolationPoint( ref_x + Point2D{ { -0.5 * ref_h, 0.5 * ref_h } },
                                                             faceStencil[vertexdof::stencilIndexFromVertex( SD::VERTEX_NW )] );
                     horiInterpolator.addInterpolationPoint( ref_x + Point2D{ { -0.5 * ref_h, 0.0 } },
                                                             faceStencil[vertexdof::stencilIndexFromVertex( SD::VERTEX_W )] );
                  }

                  x += d0;
               }

               vertInterpolator.addInterpolationPoint( ref_x + Point2D{ { 0.0, 0.5 * ref_h } },
                                                       faceStencil[vertexdof::stencilIndexFromVertex( SD::VERTEX_N )] );

               --inner_rowsize;
            }

            centerInterpolator.interpolate( facePolynomials->getPolynomialC( polyDegree ) );
            horiInterpolator.interpolate( facePolynomials->getPolynomialW( polyDegree ) );
            vertInterpolator.interpolate( facePolynomials->getPolynomialS( polyDegree ) );
            diagInterpolator.interpolate( facePolynomials->getPolynomialSE( polyDegree ) );
         }
      }
   }

   void interpolateStencilsAsymmetric( uint_t polyDegree )
   {
      typedef stencilDirection SD;
      using namespace P1Elements;

      std::vector< real_t > faceStencil( 7 );

      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;
         form.setGeometryMap( face.getGeometryMap() );

         for ( uint_t level = minLevel_; level <= maxLevel_; ++level )
         {
            auto facePolynomials = face.getData( facePolynomialIDs_[level] );
            facePolynomials->addDegree( polyDegree );

            uint_t rowsize       = levelinfo::num_microvertices_per_edge( interpolationLevel_ );
            uint_t rowsizeFine   = levelinfo::num_microvertices_per_edge( level );
            uint_t inner_rowsize = rowsize;

            VertexInterpolator interpolatorS( polyDegree, interpolationLevel_ );
            VertexInterpolator interpolatorSE( polyDegree, interpolationLevel_ );
            VertexInterpolator interpolatorW( polyDegree, interpolationLevel_ );
            VertexInterpolator interpolatorC( polyDegree, interpolationLevel_ );
            VertexInterpolator interpolatorE( polyDegree, interpolationLevel_ );
            VertexInterpolator interpolatorNW( polyDegree, interpolationLevel_ );
            VertexInterpolator interpolatorN( polyDegree, interpolationLevel_ );

            Point3D x, x0;
            x0 = face.coords[0];

            real_t ref_H = 1.0 / ( walberla::real_c( rowsize - 1 ) );
            real_t ref_h = 1.0 / ( walberla::real_c( rowsizeFine - 1 ) );

            Point3D d0 = ref_H * ( face.coords[1] - face.coords[0] );
            Point3D d2 = ref_H * ( face.coords[2] - face.coords[0] );

            // fine directions
            Point3D d0f = ref_h * ( face.coords[1] - face.coords[0] );
            Point3D d2f = ref_h * ( face.coords[2] - face.coords[0] );

            Point2D ref_x;

            Point3D dirS  = -1.0 * d2f;
            Point3D dirSE = d0f - 1.0 * d2f;
            Point3D dirE  = d0f;
            Point3D dirW  = -1.0 * d0f;
            Point3D dirNW = -1.0 * d0f + d2f;
            Point3D dirN  = d2f;

            for ( uint_t j = 1; j < rowsize - 2; ++j )
            {
               ref_x[1] = j * ref_H;

               x = x0;
               x += walberla::real_c( j ) * d2 + d0;

               uint_t i;
               for ( i = 1; i < inner_rowsize - 2; ++i )
               {
                  ref_x[0] = i * ref_H;

                  std::fill( faceStencil.begin(), faceStencil.end(), 0.0 );

                  vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                      form, { x, x + dirW, x + dirS }, P1Elements::P1Elements2D::elementSW, faceStencil.data() );
                  vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                      form, { x, x + dirS, x + dirSE }, P1Elements::P1Elements2D::elementS, faceStencil.data() );
                  vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                      form, { x, x + dirSE, x + dirE }, P1Elements::P1Elements2D::elementSE, faceStencil.data() );
                  vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                      form, { x, x + dirE, x + dirN }, P1Elements::P1Elements2D::elementNE, faceStencil.data() );
                  vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                      form, { x, x + dirN, x + dirNW }, P1Elements::P1Elements2D::elementN, faceStencil.data() );
                  vertexdof::variablestencil::assembleLocalStencil< P1Form >(
                      form, { x, x + dirNW, x + dirW }, P1Elements::P1Elements2D::elementNW, faceStencil.data() );

                  //if (i == 1 && j == 1) {
                  //   PointND<real_t, 7> test(faceStencil.data());
                  //   WALBERLA_LOG_INFO("stencil = " << test);
                  //}

                  interpolatorS.addInterpolationPoint( ref_x, faceStencil[vertexdof::stencilIndexFromVertex( SD::VERTEX_S )] );
                  interpolatorSE.addInterpolationPoint( ref_x, faceStencil[vertexdof::stencilIndexFromVertex( SD::VERTEX_SE )] );
                  interpolatorW.addInterpolationPoint( ref_x, faceStencil[vertexdof::stencilIndexFromVertex( SD::VERTEX_W )] );
                  interpolatorC.addInterpolationPoint( ref_x, faceStencil[vertexdof::stencilIndexFromVertex( SD::VERTEX_C )] );
                  interpolatorE.addInterpolationPoint( ref_x, faceStencil[vertexdof::stencilIndexFromVertex( SD::VERTEX_E )] );
                  interpolatorNW.addInterpolationPoint( ref_x, faceStencil[vertexdof::stencilIndexFromVertex( SD::VERTEX_NW )] );
                  interpolatorN.addInterpolationPoint( ref_x, faceStencil[vertexdof::stencilIndexFromVertex( SD::VERTEX_N )] );

                  x += d0;
               }
               --inner_rowsize;
            }

            interpolatorS.interpolate( facePolynomials->getPolynomialS( polyDegree ) );
            interpolatorSE.interpolate( facePolynomials->getPolynomialSE( polyDegree ) );
            interpolatorW.interpolate( facePolynomials->getPolynomialW( polyDegree ) );
            interpolatorC.interpolate( facePolynomials->getPolynomialC( polyDegree ) );
            interpolatorE.interpolate( facePolynomials->getPolynomialE( polyDegree ) );
            interpolatorNW.interpolate( facePolynomials->getPolynomialNW( polyDegree ) );
            interpolatorN.interpolate( facePolynomials->getPolynomialN( polyDegree ) );
         }
      }
   }

   void useDegree( uint_t degree ) { polyDegree_ = degree; }

   void apply( const P1Function< real_t >& src,
               const P1Function< real_t >& dst,
               const size_t                level,
               DoFType                     flag,
               UpdateType                  updateType = Replace ) const override final
   {
      WALBERLA_ASSERT_NOT_IDENTICAL( std::addressof( src ), std::addressof( dst ) );

      checkForMissingPolynomial( level, polyDegree_ );

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
            if ( OprType == OperatorType::ODD )
            {
               vertexdof::macroface::applyPolynomialOdd< real_t >( polyDegree_,
                                                                   level,
                                                                   face,
                                                                   facePolynomialIDs_.at( level ),
                                                                   src.getFaceDataID(),
                                                                   dst.getFaceDataID(),
                                                                   updateType );
            }
            else
            {
               vertexdof::macroface::applyPolynomial< real_t, OprType >( polyDegree_,
                                                                         level,
                                                                         face,
                                                                         facePolynomialIDs_.at( level ),
                                                                         src.getFaceDataID(),
                                                                         dst.getFaceDataID(),
                                                                         updateType );
            }
         }
      }
   }

   void smooth_gs( const P1Function< real_t >& dst, const P1Function< real_t >& rhs, size_t level, DoFType flag ) const override
   {
      checkForMissingPolynomial( level, polyDegree_ );

      // start pulling vertex halos
      // dst.startCommunication< Edge, Vertex >( level );

      // start pulling edge halos
      // dst.startCommunication< Face, Edge >( level );

      // end pulling vertex halos
      // dst.endCommunication< Edge, Vertex >( level );

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

      // dst.startCommunication< Vertex, Edge >( level );

      // end pulling edge halos
      // dst.endCommunication< Face, Edge >( level );

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

      // dst.endCommunication< Vertex, Edge >( level );

      // dst.startCommunication< Edge, Face >( level );

      dst.communicate< Edge, Face >( level );

      for ( auto& it : storage_->getFaces() )
      {
         Face& face = *it.second;

         const DoFType faceBC = dst.getBoundaryCondition().getBoundaryType( face.getMeshBoundaryFlag() );
         if ( testFlag( faceBC, flag ) )
         {
            switch ( OprType )
            {
            case OperatorType::MASS:
               WALBERLA_ABORT( "Not implemented" );
               break;
            case OperatorType::EVEN:
               vertexdof::macroface::smooth_gs_polynomial_even< real_t >(
                   polyDegree_, level, face, facePolynomialIDs_.at( level ), dst.getFaceDataID(), rhs.getFaceDataID() );
               break;
            case OperatorType::ODD:
               WALBERLA_ABORT( "Not implemented" );
               break;
            }
         }
      }

      // dst.endCommunication< Edge, Face >( level );
   }

   void smooth_jac( const P1Function< real_t >& dst,
                    const P1Function< real_t >& rhs,
                    const P1Function< real_t >& tmp,
                    size_t                      level,
                    DoFType                     flag ) const override
   {
      checkForMissingPolynomial( level, polyDegree_ );

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
      checkForMissingPolynomial( level, polyDegree_ );

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
   std::map< uint_t, PrimitiveDataID< FaceP1PolynomialMemory, Face > > facePolynomialIDs_;

 private:
   void checkForMissingPolynomial( uint_t level, uint_t degree ) const
   {
      WALBERLA_ASSERT( facePolynomialIDs_.count( level ) > 0, "Polynomial for level " << level << " has not been interpolated" );
   }

   uint_t polyDegree_;
   uint_t interpolationLevel_;
   P1Form form;
};

// typedef P1PolynomialBlendingOperator< P1Form_laplace, OperatorType::EVEN > P1PolynomialBlendingLaplaceOperator;
// typedef P1PolynomialBlendingOperator< P1Form_mass, OperatorType::MASS >    P1PolynomialBlendingMassOperator;

typedef P1PolynomialBlendingOperator< forms::p1_diffusion_blending_q1, OperatorType::EVEN > P1PolynomialBlendingLaplaceOperator;
typedef P1PolynomialBlendingOperator< forms::p1_mass_blending_q4, OperatorType::MASS >      P1PolynomialBlendingMassOperator;

// typedef P1PolynomialBlendingOperator< P1Form_epsilon_11, OperatorType::EVEN > P1PolynomialBlendingEpsilonOperator_11;
// typedef P1PolynomialBlendingOperator< P1Form_epsilon_12, OperatorType::EVEN > P1PolynomialBlendingEpsilonOperator_12;
// typedef P1PolynomialBlendingOperator< P1Form_epsilon_21, OperatorType::EVEN > P1PolynomialBlendingEpsilonOperator_21;
// typedef P1PolynomialBlendingOperator< P1Form_epsilon_22, OperatorType::EVEN > P1PolynomialBlendingEpsilonOperator_22;

typedef P1PolynomialBlendingOperator< forms::p1_epsiloncc_0_0_blending_q2, OperatorType::EVEN > P1PolynomialBlendingEpsilonOperator_11;
typedef P1PolynomialBlendingOperator< forms::p1_epsiloncc_0_1_blending_q2, OperatorType::EVEN > P1PolynomialBlendingEpsilonOperator_12;
typedef P1PolynomialBlendingOperator< forms::p1_epsiloncc_1_0_blending_q2, OperatorType::EVEN > P1PolynomialBlendingEpsilonOperator_21;
typedef P1PolynomialBlendingOperator< forms::p1_epsiloncc_1_1_blending_q2, OperatorType::EVEN > P1PolynomialBlendingEpsilonOperator_22;

typedef P1PolynomialBlendingOperator< forms::p1_divt_0_blending_q1 , OperatorType::ODD > P1PolynomialBlendingDivTOperator_1;
typedef P1PolynomialBlendingOperator< forms::p1_divt_1_blending_q1, OperatorType::ODD > P1PolynomialBlendingDivTOperator_2;

typedef P1PolynomialBlendingOperator< forms::p1_div_0_blending_q1, OperatorType::ODD > P1PolynomialBlendingDivOperator_1;
typedef P1PolynomialBlendingOperator< forms::p1_div_1_blending_q1, OperatorType::ODD > P1PolynomialBlendingDivOperator_2;

typedef P1PolynomialBlendingOperator< P1Form_pspg, OperatorType::EVEN > P1PolynomialBlendingPSPGOperator;

} // namespace hyteg
