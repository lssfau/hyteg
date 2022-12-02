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
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

#include "Eigen/Eigen"

namespace hyteg {
namespace dg {

class DGForm
{
 public:
   /// \brief Returns true if only volume integrals are non-zero.
   ///        The integration of facet integrals can then be skipped altogether.
   [[nodiscard]] virtual bool onlyVolumeIntegrals() const { return false; }

   /// \brief Integrates the volume contribution on a triangle element.
   ///
   /// \param dim         2 for 2D volumes, 3 for 3D volumes
   /// \param coords      coordinates of the triangle element vertices
   /// \param trialBasis  trial function basis - determines the number of columns of the element matrix
   /// \param testBasis   test function basis - determines the number of rows of the element matrix
   /// \param trialDegree polynomial degree of the trial function on this element
   /// \param testDegree  polynomial degree of the test function on this element
   /// \param elMat       computed local element matrix
   void integrateVolume( int                                                      dim,
                         const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                         const DGBasisInfo&                                       trialBasis,
                         const DGBasisInfo&                                       testBasis,
                         int                                                      trialDegree,
                         int                                                      testDegree,
                         Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      if ( dim == 2 )
      {
         integrateVolume2D( coords, trialBasis, testBasis, trialDegree, testDegree, elMat );
      }
      else if ( dim == 3 )
      {
         integrateVolume3D( coords, trialBasis, testBasis, trialDegree, testDegree, elMat );
      }
      else
      {
         WALBERLA_ABORT( "Only supporting 2D and 3D." );
      }
   };

   /// \brief Integrates the inner part of the facet contribution on a triangle element.
   ///
   /// \param dim                2 for 2D volumes, 3 for 3D volumes
   /// \param coordsElement      coordinates of the triangle element vertices
   /// \param coordsFacet        coordinates of the interface edge
   /// \param oppositeVertex     coordinates of the vertex that is part of the inner element and opposite of the interface
   /// \param outwardNormal      (normalized) normal vector to the interface in direction of the outer element
   /// \param trialBasis         trial function basis - determines the number of columns of the element matrix
   /// \param testBasis          test function basis - determines the number of rows of the element matrix
   /// \param trialDegree        polynomial degree of the trial function on this element
   /// \param testDegree         polynomial degree of the test function on this element
   /// \param elMat              computed local element matrix
   void integrateFacetInner( int                                                      dim,
                             const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                             const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                             const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                             const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                             const DGBasisInfo&                                       trialBasis,
                             const DGBasisInfo&                                       testBasis,
                             int                                                      trialDegree,
                             int                                                      testDegree,
                             Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      if ( dim == 2 )
      {
         integrateFacetInner2D(
             coordsElement, coordsFacet, oppositeVertex, outwardNormal, trialBasis, testBasis, trialDegree, testDegree, elMat );
      }
      else if ( dim == 3 )
      {
         integrateFacetInner3D(
             coordsElement, coordsFacet, oppositeVertex, outwardNormal, trialBasis, testBasis, trialDegree, testDegree, elMat );
      }
      else
      {
         WALBERLA_ABORT( "Only supporting 2D and 3D." );
      }
   }

   /// \brief Integrates the facet contributions from both sides on a triangle element - writing to the inner DoFs.
   ///
   /// \param dim                        2 for 2D volumes, 3 for 3D volumes
   /// \param coordsElementInner         coordinates of the triangle element vertices of the inner (dst) element
   /// \param coordsElementOuter         coordinates of the triangle element vertices of the outer (src) element
   /// \param coordsFacet                coordinates of the interface edge
   /// \param oppositeVertexInnerElement coordinates of the vertex that is part of the inner element and opposite of the interface
   /// \param oppositeVertexOuterElement coordinates of the vertex that is part of the outer element and opposite of the interface
   /// \param outwardNormal              (normalized) normal vector to the interface in direction of the outer element
   /// \param trialBasis                 trial function basis - determines the number of columns of the element matrix (outer element)
   /// \param testBasis                  test function basis - determines the number of rows of the element matrix (inner element)
   /// \param trialDegree                polynomial degree of the trial function on this element
   /// \param testDegree                 polynomial degree of the test function on this element
   /// \param elMat                      computed local element matrix
   void integrateFacetCoupling( int                                                      dim,
                                const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                const DGBasisInfo&                                       trialBasis,
                                const DGBasisInfo&                                       testBasis,
                                int                                                      trialDegree,
                                int                                                      testDegree,
                                Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      if ( dim == 2 )
      {
         integrateFacetCoupling2D( coordsElementInner,
                                   coordsElementOuter,
                                   coordsFacet,
                                   oppositeVertexInnerElement,
                                   oppositeVertexOuterElement,
                                   outwardNormal,
                                   trialBasis,
                                   testBasis,
                                   trialDegree,
                                   testDegree,
                                   elMat );
      }
      else if ( dim == 3 )
      {
         integrateFacetCoupling3D( coordsElementInner,
                                   coordsElementOuter,
                                   coordsFacet,
                                   oppositeVertexInnerElement,
                                   oppositeVertexOuterElement,
                                   outwardNormal,
                                   trialBasis,
                                   testBasis,
                                   trialDegree,
                                   testDegree,
                                   elMat );
      }
      else
      {
         WALBERLA_ABORT( "Only supporting 2D and 3D." );
      }
   }

   /// \brief Integrates the facet contributions at Dirichlet boundaries - writing to the inner DoFs.
   ///
   /// \param dim                     2 for 2D volumes, 3 for 3D volumes
   /// \param coordsElement           coordinates of the triangle element vertices of the inner (dst) element
   /// \param coordsFacet             coordinates of the interface edge
   /// \param oppositeVertex          coordinates of the vertex that is part of the inner element and opposite of the interface
   /// \param outwardNormal           (normalized) normal vector to the interface in direction of the outer element
   /// \param trialBasis              trial function basis - determines the number of columns of the element matrix (outer element)
   /// \param testBasis               test function basis - determines the number of rows of the element matrix (inner element)
   /// \param trialDegree             polynomial degree of the trial function on this element
   /// \param testDegree              polynomial degree of the test function on this element
   /// \param elMat                   computed local element matrix
   void integrateFacetDirichletBoundary( int                                                      dim,
                                         const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                         const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                         const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                         const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                         const DGBasisInfo&                                       trialBasis,
                                         const DGBasisInfo&                                       testBasis,
                                         int                                                      trialDegree,
                                         int                                                      testDegree,
                                         Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      if ( dim == 2 )
      {
         integrateFacetDirichletBoundary2D(
             coordsElement, coordsFacet, oppositeVertex, outwardNormal, trialBasis, testBasis, trialDegree, testDegree, elMat );
      }
      else if ( dim == 3 )
      {
         integrateFacetDirichletBoundary3D(
             coordsElement, coordsFacet, oppositeVertex, outwardNormal, trialBasis, testBasis, trialDegree, testDegree, elMat );
      }
      else
      {
         WALBERLA_ABORT( "Only supporting 2D and 3D." );
      }
   }

   /// \brief Integrates the facet contributions at Dirichlet boundaries for the right-hand side.
   ///
   /// \param dim                     2 for 2D volumes, 3 for 3D volumes
   /// \param coordsElement           coordinates of the triangle element vertices of the inner (dst) element
   /// \param coordsFacet             coordinates of the interface edge
   /// \param oppositeVertex          coordinates of the vertex that is part of the inner element and opposite of the interface
   /// \param outwardNormal           (normalized) normal vector to the interface in direction of the outer element
   /// \param asis                    function basis - determines the number of rows of the element matrix
   /// \param degree                  polynomial degree of the function on this element
   /// \param elMat                   computed local element matrix - returned as a column vector
   void integrateRHSDirichletBoundary( int                                                      dim,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       basis,
                                       int                                                      degree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const
   {
      if ( dim == 2 )
      {
         integrateRHSDirichletBoundary2D( coordsElement, coordsFacet, oppositeVertex, outwardNormal, basis, degree, elMat );
      }
      else if ( dim == 3 )
      {
         integrateRHSDirichletBoundary3D( coordsElement, coordsFacet, oppositeVertex, outwardNormal, basis, degree, elMat );
      }
      else
      {
         WALBERLA_ABORT( "Only supporting 2D and 3D." );
      }
   }

 protected:
   virtual void integrateVolume2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                                   const DGBasisInfo&                                       trialBasis,
                                   const DGBasisInfo&                                       testBasis,
                                   int                                                      trialDegree,
                                   int                                                      testDegree,
                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const = 0;

   virtual void integrateVolume3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coords,
                                   const DGBasisInfo&                                       trialBasis,
                                   const DGBasisInfo&                                       testBasis,
                                   int                                                      trialDegree,
                                   int                                                      testDegree,
                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const = 0;

   virtual void integrateFacetInner2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const = 0;

   virtual void integrateFacetInner3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                       const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                       const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                       const DGBasisInfo&                                       trialBasis,
                                       const DGBasisInfo&                                       testBasis,
                                       int                                                      trialDegree,
                                       int                                                      testDegree,
                                       Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const = 0;

   virtual void integrateFacetCoupling2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const = 0;

   virtual void integrateFacetCoupling3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementInner,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElementOuter,
                                          const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexInnerElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertexOuterElement,
                                          const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                          const DGBasisInfo&                                       trialBasis,
                                          const DGBasisInfo&                                       testBasis,
                                          int                                                      trialDegree,
                                          int                                                      testDegree,
                                          Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const = 0;

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const = 0;

   virtual void integrateFacetDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                   const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                   const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                   const DGBasisInfo&                                       trialBasis,
                                                   const DGBasisInfo&                                       testBasis,
                                                   int                                                      trialDegree,
                                                   int                                                      testDegree,
                                                   Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const = 0;

   virtual void integrateRHSDirichletBoundary2D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const = 0;

   virtual void integrateRHSDirichletBoundary3D( const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsElement,
                                                 const std::vector< Eigen::Matrix< real_t, 3, 1 > >&      coordsFacet,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     oppositeVertex,
                                                 const Eigen::Matrix< real_t, 3, 1 >&                     outwardNormal,
                                                 const DGBasisInfo&                                       basis,
                                                 int                                                      degree,
                                                 Eigen::Matrix< real_t, Eigen::Dynamic, Eigen::Dynamic >& elMat ) const = 0;
};

} // namespace dg
} // namespace hyteg