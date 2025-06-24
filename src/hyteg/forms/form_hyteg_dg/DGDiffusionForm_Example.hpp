/*
* Copyright (c) 2017-2025 Nils Kohl, Marcus Mohr.
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

#include "core/DataTypes.h"

#include "hyteg/dgfunctionspace/DGBasisInfo.hpp"
#include "hyteg/forms/DGForm.hpp"
#include "hyteg/forms/form_hyteg_dg/DGForm2D.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

#include "Eigen/Eigen"

namespace hyteg {
namespace dg {

using walberla::real_c;

/// \brief Reference implementation of a DG form for the Laplacian.
///
/// Refer to
///     Béatrice Rivière (2008).
///     Discontinuous Galerkin Methods for Solving Elliptic and Parabolic Equations: Theory and Implementation
/// for details.
///
/// Note: if uncertain, try setting beta_0 to 1 for 2D and 0.5 for 3D applications.
class DGDiffusionForm_Example : public DGForm
{
 public:
   DGDiffusionForm_Example( real_t _beta_0 )
   : callback_Scalar_Variable_Coefficient_2D_g( []( const Point3D& ) { return real_c( 0 ); } )
   , callback_Scalar_Variable_Coefficient_3D_g( []( const Point3D& ) { return real_c( 0 ); } )
   , sigma_0( 6 )
   , beta_0( _beta_0 )
   {}

   DGDiffusionForm_Example( real_t                                           _beta_0,
                            const std::function< real_t( const Point3D& ) >& _callback_Scalar_Variable_Coefficient_2D_g,
                            const std::function< real_t( const Point3D& ) >& _callback_Scalar_Variable_Coefficient_3D_g )
   : callback_Scalar_Variable_Coefficient_2D_g( _callback_Scalar_Variable_Coefficient_2D_g )
   , callback_Scalar_Variable_Coefficient_3D_g( _callback_Scalar_Variable_Coefficient_3D_g )
   , sigma_0( 6 )
   , beta_0( _beta_0 )
   {}

 protected:
   void integrateVolume2D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override;

   void integrateVolume3D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&            trialBasis,
                           const DGBasisInfo&            testBasis,
                           int                           trialDegree,
                           int                           testDegree,
                           MatrixXr&                     elMat ) const override;

   void integrateFacetInner2D( const std::vector< Point3D >& coordsElement,
                               const std::vector< Point3D >& coordsFacet,
                               const Point3D&                oppositeVertex,
                               const Point3D&                outwardNormal,
                               const DGBasisInfo&            trialBasis,
                               const DGBasisInfo&            testBasis,
                               int                           trialDegree,
                               int                           testDegree,
                               MatrixXr&                     elMat ) const override;

   void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                               const std::vector< Point3D >& coordsFacet,
                               const Point3D&                oppositeVertex,
                               const Point3D&                outwardNormal,
                               const DGBasisInfo&            trialBasis,
                               const DGBasisInfo&            testBasis,
                               int                           trialDegree,
                               int                           testDegree,
                               MatrixXr&                     elMat ) const override;

   void integrateFacetCoupling2D( const std::vector< Point3D >& coordsElementInner,
                                  const std::vector< Point3D >& coordsElementOuter,
                                  const std::vector< Point3D >& coordsFacet,
                                  const Point3D&                oppositeVertexInnerElement,
                                  const Point3D&                oppositeVertexOuterElement,
                                  const Point3D&                outwardNormal,
                                  const DGBasisInfo&            trialBasis,
                                  const DGBasisInfo&            testBasis,
                                  int                           trialDegree,
                                  int                           testDegree,
                                  MatrixXr&                     elMat ) const override;

   void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                  const std::vector< Point3D >& coordsElementOuter,
                                  const std::vector< Point3D >& coordsFacet,
                                  const Point3D&                oppositeVertexInnerElement,
                                  const Point3D&                oppositeVertexOuterElement,
                                  const Point3D&                outwardNormal,
                                  const DGBasisInfo&            trialBasis,
                                  const DGBasisInfo&            testBasis,
                                  int                           trialDegree,
                                  int                           testDegree,
                                  MatrixXr&                     elMat ) const override;

   void integrateFacetDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                           const std::vector< Point3D >& coordsFacet,
                                           const Point3D&                oppositeVertex,
                                           const Point3D&                outwardNormal,
                                           const DGBasisInfo&            trialBasis,
                                           const DGBasisInfo&            testBasis,
                                           int                           trialDegree,
                                           int                           testDegree,
                                           MatrixXr&                     elMat ) const override;

   void integrateFacetDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                           const std::vector< Point3D >& coordsFacet,
                                           const Point3D&                oppositeVertex,
                                           const Point3D&                outwardNormal,
                                           const DGBasisInfo&            trialBasis,
                                           const DGBasisInfo&            testBasis,
                                           int                           trialDegree,
                                           int                           testDegree,
                                           MatrixXr&                     elMat ) const override;

   void integrateRHSDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override;

   void integrateRHSDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                         const std::vector< Point3D >& coordsFacet,
                                         const Point3D&                oppositeVertex,
                                         const Point3D&                outwardNormal,
                                         const DGBasisInfo&            basis,
                                         int                           degree,
                                         MatrixXr&                     elMat ) const override;

 private:
   void Scalar_Variable_Coefficient_2D_g( real_t in_0, real_t in_1, real_t* out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_g( Point3D( in_0, in_1, 0 ) );
   }

   void Scalar_Variable_Coefficient_3D_g( real_t in_0, real_t in_1, real_t in_2, real_t* out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_3D_g( Point3D( in_0, in_1, in_2 ) );
   }

   std::function< real_t( const Point3D& ) > callback_Scalar_Variable_Coefficient_2D_g;
   std::function< real_t( const Point3D& ) > callback_Scalar_Variable_Coefficient_3D_g;

   real_t sigma_0;
   real_t beta_0;
};

} // namespace dg
} // namespace hyteg
