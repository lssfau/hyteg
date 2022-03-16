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

#include <hyteg/operators/Operator.hpp>

#include "hyteg/types/pointnd.hpp"

#include <hyteg/p2functionspace/polynomial/P2PolynomialDataHandling.hpp>

#include "variablestencil/P2VariableStencilCommon.hpp"

#include <hyteg/forms/form_hyteg_manual/P2FormLaplace.hpp>
#include <hyteg/forms/form_hyteg_manual/P2FormDivKGrad.hpp>
#include <hyteg/forms/form_hyteg_generated/p2/p2_mass_blending_q4.hpp>

#include <hyteg/p2functionspace/polynomial/StencilInterpolator.hpp>

#include <hyteg/p1functionspace/VertexDoFFunction.hpp>
#include <hyteg/p2functionspace/P2Function.hpp>
#include <hyteg/p2functionspace/polynomial/P2MacroFacePolynomial.hpp>

#include "hyteg/communication/Syncing.hpp"


namespace hyteg {

template <class P2Form, OperatorType OprType>
class P2SurrogateOperator : public Operator<P2Function<real_t>, P2Function<real_t>>,
                            public GSSmoothable< P2Function< real_t > >,
                            public ConstantJacobiSmoothable< P1Function< real_t > >
{
public:
P2SurrogateOperator(const std::shared_ptr<PrimitiveStorage>& storage,
                                uint_t minLevel, uint_t maxLevel, uint_t interpolationLevel)
      : Operator(storage, minLevel, maxLevel)
      , interpolationLevel_(interpolationLevel)
   {
      for (uint_t level = minLevel_; level <= maxLevel_; ++level)
      {
         PrimitiveDataID<P2::FacePolynomialMemory, Face> id;
         auto dataHandling = std::make_shared<FaceP2PolynomialMemoryDataHandling>();
         storage_->addFaceData(id, dataHandling, "P2OperatorFacePolynomial");
         polynomialIDs_[level] = id;
      }
   }

   P2SurrogateOperator(const std::shared_ptr<PrimitiveStorage>& storage,
                                uint_t minLevel, uint_t maxLevel, uint_t interpolationLevel, uint_t polyDegree)
      : P2SurrogateOperator(storage, minLevel, maxLevel, interpolationLevel)
   {
      interpolateStencils(polyDegree);
      useDegree(polyDegree);
   }

   ~P2SurrogateOperator() {}

   void interpolateStencils(uint_t polyDegree)
   {
      real_t  H = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(interpolationLevel_) - 1));

      // stencil entries
      std::array<real_t, P2::NumStencilentries2D::VtV> vertexToVertexStencil;
      std::array<real_t, P2::NumStencilentries2D::EtV> edgeToVertexStencil;
      std::array<real_t, P2::NumStencilentries2D::VtE> vertexToEdgeStencil;
      std::array<real_t, P2::NumStencilentries2D::EtE> edgeToEdgeStencil;

      // we only use polynomials for face stencils
      for (auto& itF : storage_->getFaces())
      {
         Face& face = *itF.second;
         form_.setGeometryMap(face.getGeometryMap());

         Point3D x0(face.coords[0]), x;

         Point3D D0 = face.coords[1] - face.coords[0];
         Point3D D2 = face.coords[2] - face.coords[0];

         for (uint_t level = minLevel_; level <= maxLevel_; ++level)
         {
            P2::StencilInterpolator<LSQPType::VERTEX, P2::NumStencilentries2D::VtV> VtVInterpolator(polyDegree, interpolationLevel_);
            P2::StencilInterpolator<LSQPType::VERTEX, P2::NumStencilentries2D::EtV> EtVInterpolator(polyDegree, interpolationLevel_);
            P2::StencilInterpolator<LSQPType::EDGE_ALL, P2::NumStencilentries2D::VtE> VtEInterpolator(polyDegree, interpolationLevel_);
            P2::StencilInterpolator<LSQPType::EDGE_ALL, P2::NumStencilentries2D::EtE> EtEInterpolator(polyDegree, interpolationLevel_);

            // directions (size of microfaces based on current level)
            real_t  h = 1.0 / (walberla::real_c(levelinfo::num_microvertices_per_edge(level) - 1));
            const Point3D dirS  = h * (-D2);
            const Point3D dirSE = h * (D0 - D2);
            const Point3D dirE  = h * (D0);
            const Point3D dirW  = -dirE;
            const Point3D dirNW = -dirSE;
            const Point3D dirN  = -dirS;
            const Point3D dirNE = dirN + dirE;

            // loop over all DOFs (number and position of microfaces based of interpolationLevel)
            for (const auto& it : hyteg::edgedof::macroface::Iterator(interpolationLevel_, 0))
            {
               // position of vertex on macroface
               real_t r = walberla::real_c(it.row()) * H;
               real_t c = walberla::real_c(it.col()) * H;
               x = x0 + (r * D2 + c * D0);
               // corresponding point on reference element
               Point2D xi({c, r});

               P2::variablestencil::macroface::assembleStencil(form_, x, dirS, dirSE, dirE, dirN, dirNW, dirW, dirNE, vertexToVertexStencil, edgeToVertexStencil, vertexToEdgeStencil, edgeToEdgeStencil);

               // add interpolation points
               // vertex DoF
               if (!vertexdof::macroface::isVertexOnBoundary(interpolationLevel_, it))
               {
                  for (uint_t i = 0; i < P2::NumStencilentries2D::VtV; ++i)
                  {
                     VtVInterpolator[i].addInterpolationPoint(xi, vertexToVertexStencil[i]);
                  }

                  for (uint_t i = 0; i < P2::NumStencilentries2D::EtV; ++i)
                  {
                     EtVInterpolator[i].addInterpolationPoint(xi, edgeToVertexStencil[i]);
                  }
               }

               // horizontal edge DoF
               if (!edgedof::macroface::isHorizontalEdgeOnBoundary(interpolationLevel_, it))
               {
                  for (const auto& dir : vertexdof::macroface::neighborsFromHorizontalEdge)
                  {
                     const uint_t i = vertexdof::stencilIndexFromHorizontalEdge(dir);
                     VtEInterpolator[i].addInterpolationPoint(xi, vertexToEdgeStencil[i]);
                  }

                  for (const auto& dir : edgedof::macroface::neighborsFromHorizontalEdge)
                  {
                     const uint_t i = edgedof::stencilIndexFromHorizontalEdge(dir);
                     EtEInterpolator[i].addInterpolationPoint(xi, edgeToEdgeStencil[i]);
                  }
               }

               // vertical edge DoF
               if (!edgedof::macroface::isVerticalEdgeOnBoundary(interpolationLevel_, it))
               {
                  for (const auto& dir : vertexdof::macroface::neighborsFromVerticalEdge)
                  {
                     const uint_t i = vertexdof::stencilIndexFromVerticalEdge(dir);
                     VtEInterpolator[i].addInterpolationPoint(xi, vertexToEdgeStencil[i]);
                  }

                  for (const auto& dir : edgedof::macroface::neighborsFromVerticalEdge)
                  {
                     const uint_t i = edgedof::stencilIndexFromVerticalEdge(dir);
                     EtEInterpolator[i].addInterpolationPoint(xi, edgeToEdgeStencil[i]);
                  }
               }

               // diagonal edge DoF
               if (!edgedof::macroface::isDiagonalEdgeOnBoundary(interpolationLevel_, it))
               {
                  for (const auto& dir : vertexdof::macroface::neighborsFromDiagonalEdge)
                  {
                     const uint_t i = vertexdof::stencilIndexFromDiagonalEdge(dir);
                     VtEInterpolator[i].addInterpolationPoint(xi, vertexToEdgeStencil[i]);
                  }

                  for (const auto& dir : edgedof::macroface::neighborsFromDiagonalEdge)
                  {
                     const uint_t i = edgedof::stencilIndexFromDiagonalEdge(dir);
                     EtEInterpolator[i].addInterpolationPoint(xi, edgeToEdgeStencil[i]);
                  }
               }
            }

            auto id = face.getData(polynomialIDs_[level]);

            auto& poly = id->addDegree(polyDegree);

            VtVInterpolator.interpolate(poly.VtV);
            EtVInterpolator.interpolate(poly.EtV);
            VtEInterpolator.interpolate(poly.VtE);
            EtEInterpolator.interpolate(poly.EtE);
         }
      }
   }

   void useDegree(uint_t degree) {polyDegree_ = degree;}

   void apply(const P2Function<real_t>& src, const P2Function<real_t>& dst,
              const size_t level, DoFType flag, UpdateType updateType = Replace) const
   {
      WALBERLA_ASSERT_NOT_IDENTICAL(std::addressof(src), std::addressof(dst));

      checkForMissingPolynomial(level);

      communication::syncP2FunctionBetweenPrimitives(src, level);

      const vertexdof::VertexDoFFunction<real_t>&  srcVertexDoF   = src.getVertexDoFFunction();
      const EdgeDoFFunction<real_t>&               srcEdgeDoF     = src.getEdgeDoFFunction();
      const vertexdof::VertexDoFFunction<real_t>&  dstVertexDoF   = dst.getVertexDoFFunction();
      const EdgeDoFFunction<real_t>&               dstEdgeDoF     = dst.getEdgeDoFFunction();

      for (auto& it : storage_->getVertices())
      {
         hyteg::Vertex& vertex = *it.second;

         const DoFType vtxFlag = dst.getBoundaryCondition().getBoundaryType(vertex.getMeshBoundaryFlag());

         if (testFlag(vtxFlag, flag))
         {
            P2::variablestencil::macrovertex::applyVariableStencil<P2Form>(level, vertex, storage_, srcVertexDoF.getVertexDataID(), srcEdgeDoF.getVertexDataID(), dstVertexDoF.getVertexDataID(), updateType);
         }
      }

      for (auto& it : storage_->getEdges())
      {
         hyteg::Edge& edge = *it.second;

         const DoFType edgeFlag = dst.getBoundaryCondition().getBoundaryType(edge.getMeshBoundaryFlag());

         if (testFlag(edgeFlag, flag))
         {
            P2::variablestencil::macroedge::applyVariableStencil<P2Form>(level, edge, storage_, srcVertexDoF.getEdgeDataID(), srcEdgeDoF.getEdgeDataID(), dstVertexDoF.getEdgeDataID(), dstEdgeDoF.getEdgeDataID(), updateType);
         }
      }

      for (auto& it : storage_->getFaces())
      {
         hyteg::Face& face = *it.second;

         const DoFType faceFlag = dst.getBoundaryCondition().getBoundaryType(face.getMeshBoundaryFlag());

         if (testFlag(faceFlag, flag))
         {
            P2::variablestencil::macroface::applyPolynomial(
               polyDegree_, polynomialIDs_.at(level), level, face,
               srcVertexDoF.getFaceDataID(), srcEdgeDoF.getFaceDataID(),
               dstVertexDoF.getFaceDataID(), dstEdgeDoF.getFaceDataID(), updateType);
         }
      }
   }

   void smooth_gs(const P2Function<real_t>& dst, const P2Function<real_t>& rhs,
                  const size_t level, DoFType flag) const override
   {
      checkForMissingPolynomial(level);

      communication::syncP2FunctionBetweenPrimitives(dst, level);

      const vertexdof::VertexDoFFunction<real_t>&  dstVertexDoF   = dst.getVertexDoFFunction();
      const EdgeDoFFunction<real_t>&               dstEdgeDoF     = dst.getEdgeDoFFunction();
      const vertexdof::VertexDoFFunction<real_t>&  rhsVertexDoF   = rhs.getVertexDoFFunction();
      const EdgeDoFFunction<real_t>&               rhsEdgeDoF     = rhs.getEdgeDoFFunction();

      for (auto& it : storage_->getVertices())
      {
         hyteg::Vertex& vertex = *it.second;

         const DoFType vertexFlag = dst.getBoundaryCondition().getBoundaryType(vertex.getMeshBoundaryFlag());

         if (testFlag(vertexFlag, flag))
         {
            P2::variablestencil::macrovertex::smoothGSVariableStencil<P2Form>(level, vertex, storage_, dstVertexDoF.getVertexDataID(), dstEdgeDoF.getVertexDataID(), rhsVertexDoF.getVertexDataID());
         }
      }

      communication::syncP2FunctionBetweenPrimitives(dst, level);

      for (auto& it : storage_->getEdges())
      {
         hyteg::Edge& edge = *it.second;

         const DoFType edgeFlag = dst.getBoundaryCondition().getBoundaryType(edge.getMeshBoundaryFlag());

         if (testFlag(edgeFlag, flag))
         {
            P2::variablestencil::macroedge::smoothGSVariableStencil<P2Form>(level, edge, storage_, dstVertexDoF.getEdgeDataID(), dstEdgeDoF.getEdgeDataID(), rhsVertexDoF.getEdgeDataID(), rhsEdgeDoF.getEdgeDataID());
         }
      }

      communication::syncP2FunctionBetweenPrimitives(dst, level);

      for (auto& it : storage_->getFaces())
      {
         hyteg::Face& face = *it.second;

         const DoFType faceFlag = dst.getBoundaryCondition().getBoundaryType(face.getMeshBoundaryFlag());

         if (testFlag(faceFlag, flag))
         {
            P2::variablestencil::macroface::smoothGSPolynomial(
               polyDegree_, polynomialIDs_.at(level), level, face,
               dstVertexDoF.getFaceDataID(), dstEdgeDoF.getFaceDataID(),
               rhsVertexDoF.getFaceDataID(), rhsEdgeDoF.getFaceDataID());
         }
      }

      // communication::syncP2FunctionBetweenPrimitives(dst, level);

      if (storage_->hasGlobalCells())
      {
         WALBERLA_ABORT("P2PolynomialBleningOperator not implemented for 3D")
      }
   }

   void smooth_jac(const P1Function<real_t>& dst, const P1Function<real_t>& rhs,
                   const P1Function<real_t>& tmp, size_t level, DoFType flag) const override
   {
      WALBERLA_ABORT("To be implemented");
   }

 private:

   inline void checkForMissingPolynomial(uint_t level) const
   {
      WALBERLA_ASSERT(polynomialIDs_.count(level) > 0, "Polynomial for level " << level << " has not been interpolated");
   }

   uint_t polyDegree_;
   uint_t interpolationLevel_;
   P2Form form_;
   std::map<uint_t, PrimitiveDataID<P2::FacePolynomialMemory, Face>> polynomialIDs_;
};

typedef P2SurrogateOperator<P2Form_laplace, OperatorType::EVEN> P2SurrogateLaplaceOperator;
typedef P2SurrogateOperator<P2Form_divKgrad, OperatorType::EVEN> P2SurrogateDivKgradOperator;
typedef P2SurrogateOperator<forms::p2_mass_blending_q4, OperatorType::MASS>    P2SurrogateMassOperator;

} // namespace hyteg
