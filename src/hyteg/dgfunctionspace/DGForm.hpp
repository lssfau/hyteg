/*
* Copyright (c) 2017-2022 Nils Kohl.
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
#include "hyteg/types/matrix.hpp"
#include "hyteg/types/pointnd.hpp"

#include "Eigen/Eigen"

namespace hyteg {
namespace dg {

class DGForm
{
 public:
   /// \brief Integrates the volume contribution on a triangle element (2D).
   ///
   /// \param coords      coordinates of the triangle element vertices
   /// \param trialBasis  trial function basis - determines the number of columns of the element matrix
   /// \param testBasis   test function basis - determines the number of rows of the element matrix
   /// \param trialDegree polynomial degree of the trial function on this element
   /// \param testDegree  polynomial degree of the test function on this element
   /// \param elMat       computed local element matrix
   virtual void integrateVolume( const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >&    coords,
                                 const DGBasisInfo&                                       trialBasis,
                                 const DGBasisInfo&                                       testBasis,
                                 int                                                      trialDegree,
                                 int                                                      testDegree,
                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const = 0;

   /// \brief Integrates the inner part of the facet contribution on a triangle element (2D).
   ///
   /// \param coordsElement      coordinates of the triangle element vertices
   /// \param coordsFacet        coordinates of the interface edge
   /// \param outwardNormal      (normalized) normal vector to the interface in direction of the outer element
   /// \param trialBasis         trial function basis - determines the number of columns of the element matrix
   /// \param testBasis          test function basis - determines the number of rows of the element matrix
   /// \param trialDegree        polynomial degree of the trial function on this element
   /// \param testDegree         polynomial degree of the test function on this element
   /// \param elMat              computed local element matrix
   virtual void integrateFacetInner( const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >&    coordsElement,
                                     const std::array< Eigen::Matrix< real_t, 2, 1 >, 2 >&    coordsFacet,
                                     const Eigen::Matrix< real_t, 2, 1 >&                     outwardNormal,
                                     const DGBasisInfo&                                       trialBasis,
                                     const DGBasisInfo&                                       testBasis,
                                     int                                                      trialDegree,
                                     int                                                      testDegree,
                                     Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const = 0;

   /// \brief Integrates the facet contributions from both sides on a triangle element (2D) - writing to the inner DoFs.
   ///
   /// \param coordsElementInner      coordinates of the triangle element vertices of the inner (dst) element
   /// \param coordsElementOuter      coordinates of the triangle element vertices of the outer (src) element
   /// \param coordsFacet             coordinates of the interface edge
   /// \param outwardNormal           (normalized) normal vector to the interface in direction of the outer element
   /// \param trialBasis              trial function basis - determines the number of columns of the element matrix (outer element)
   /// \param testBasis               test function basis - determines the number of rows of the element matrix (inner element)
   /// \param trialDegree             polynomial degree of the trial function on this element
   /// \param testDegree              polynomial degree of the test function on this element
   /// \param elMat                   computed local element matrix
   virtual void integrateFacetCoupling( const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >&    coordsElementInner,
                                        const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >&    coordsElementOuter,
                                        const std::array< Eigen::Matrix< real_t, 2, 1 >, 2 >&    coordsFacet,
                                        const Eigen::Matrix< real_t, 2, 1 >&                     outwardNormal,
                                        const DGBasisInfo&                                       trialBasis,
                                        const DGBasisInfo&                                       testBasis,
                                        int                                                      trialDegree,
                                        int                                                      testDegree,
                                        Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const = 0;

   /// \brief Integrates the facet contributions at Dirichlet boundaries (2D) - writing to the inner DoFs.
   ///
   /// \param coordsElement           coordinates of the triangle element vertices of the inner (dst) element
   /// \param coordsFacet             coordinates of the interface edge
   /// \param outwardNormal           (normalized) normal vector to the interface in direction of the outer element
   /// \param trialBasis              trial function basis - determines the number of columns of the element matrix (outer element)
   /// \param testBasis               test function basis - determines the number of rows of the element matrix (inner element)
   /// \param trialDegree             polynomial degree of the trial function on this element
   /// \param testDegree              polynomial degree of the test function on this element
   /// \param elMat                   computed local element matrix
   virtual void integrateFacetDirichletBoundary( const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >&    coordsElement,
                                                 const std::array< Eigen::Matrix< real_t, 2, 1 >, 2 >&    coordsFacet,
                                                 const Eigen::Matrix< real_t, 2, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       trialBasis,
                                                 const DGBasisInfo&                                       testBasis,
                                                 int                                                      trialDegree,
                                                 int                                                      testDegree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const = 0;
};

} // namespace dg
} // namespace hyteg