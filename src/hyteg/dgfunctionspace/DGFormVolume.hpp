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
#include "hyteg/dgfunctionspace/DGForm.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

#include "Eigen/Eigen"

namespace hyteg {
namespace dg {

/// \brief Helper class to derive from if only volume integrals are non-zero.
class DGFormVolume : public DGForm
{
 public:
   virtual bool onlyVolumeIntegrals() const { return true; }

 protected:
   virtual void integrateFacetInner2D( const std::vector< Point3D >& coordsElement,
                                       const std::vector< Point3D >& coordsFacet,
                                       const Point3D&                oppositeVertex,
                                       const Point3D&                outwardNormal,
                                       const DGBasisInfo&            trialBasis,
                                       const DGBasisInfo&            testBasis,
                                       int                           trialDegree,
                                       int                           testDegree,
                                       MatrixXr&                     elMat ) const
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( trialBasis );
      WALBERLA_UNUSED( testBasis );
      WALBERLA_UNUSED( trialDegree );
      WALBERLA_UNUSED( testDegree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }

   virtual void integrateFacetInner3D( const std::vector< Point3D >& coordsElement,
                                       const std::vector< Point3D >& coordsFacet,
                                       const Point3D&                oppositeVertex,
                                       const Point3D&                outwardNormal,
                                       const DGBasisInfo&            trialBasis,
                                       const DGBasisInfo&            testBasis,
                                       int                           trialDegree,
                                       int                           testDegree,
                                       MatrixXr&                     elMat ) const
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( trialBasis );
      WALBERLA_UNUSED( testBasis );
      WALBERLA_UNUSED( trialDegree );
      WALBERLA_UNUSED( testDegree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >& coordsElementInner,
                                          const std::vector< Point3D >& coordsElementOuter,
                                          const std::vector< Point3D >& coordsFacet,
                                          const Point3D&                oppositeVertexInnerElement,
                                          const Point3D&                oppositeVertexOuterElement,
                                          const Point3D&                outwardNormal,
                                          const DGBasisInfo&            trialBasis,
                                          const DGBasisInfo&            testBasis,
                                          int                           trialDegree,
                                          int                           testDegree,
                                          MatrixXr&                     elMat ) const
   {
      WALBERLA_UNUSED( coordsElementInner );
      WALBERLA_UNUSED( coordsElementOuter );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertexInnerElement );
      WALBERLA_UNUSED( oppositeVertexOuterElement );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( trialBasis );
      WALBERLA_UNUSED( testBasis );
      WALBERLA_UNUSED( trialDegree );
      WALBERLA_UNUSED( testDegree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   };

   virtual void integrateFacetCoupling3D( const std::vector< Point3D >& coordsElementInner,
                                          const std::vector< Point3D >& coordsElementOuter,
                                          const std::vector< Point3D >& coordsFacet,
                                          const Point3D&                oppositeVertexInnerElement,
                                          const Point3D&                oppositeVertexOuterElement,
                                          const Point3D&                outwardNormal,
                                          const DGBasisInfo&            trialBasis,
                                          const DGBasisInfo&            testBasis,
                                          int                           trialDegree,
                                          int                           testDegree,
                                          MatrixXr&                     elMat ) const
   {
      WALBERLA_UNUSED( coordsElementInner );
      WALBERLA_UNUSED( coordsElementOuter );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertexInnerElement );
      WALBERLA_UNUSED( oppositeVertexOuterElement );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( trialBasis );
      WALBERLA_UNUSED( testBasis );
      WALBERLA_UNUSED( trialDegree );
      WALBERLA_UNUSED( testDegree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                                   const std::vector< Point3D >& coordsFacet,
                                                   const Point3D&                oppositeVertex,
                                                   const Point3D&                outwardNormal,
                                                   const DGBasisInfo&            trialBasis,
                                                   const DGBasisInfo&            testBasis,
                                                   int                           trialDegree,
                                                   int                           testDegree,
                                                   MatrixXr&                     elMat ) const
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( trialBasis );
      WALBERLA_UNUSED( testBasis );
      WALBERLA_UNUSED( trialDegree );
      WALBERLA_UNUSED( testDegree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }

   virtual void integrateFacetDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                                   const std::vector< Point3D >& coordsFacet,
                                                   const Point3D&                oppositeVertex,
                                                   const Point3D&                outwardNormal,
                                                   const DGBasisInfo&            trialBasis,
                                                   const DGBasisInfo&            testBasis,
                                                   int                           trialDegree,
                                                   int                           testDegree,
                                                   MatrixXr&                     elMat ) const
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( trialBasis );
      WALBERLA_UNUSED( testBasis );
      WALBERLA_UNUSED( trialDegree );
      WALBERLA_UNUSED( testDegree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }

   virtual void integrateRHSDirichletBoundary2D( const std::vector< Point3D >& coordsElement,
                                                 const std::vector< Point3D >& coordsFacet,
                                                 const Point3D&                oppositeVertex,
                                                 const Point3D&                outwardNormal,
                                                 const DGBasisInfo&            basis,
                                                 int                           degree,
                                                 MatrixXr&                     elMat ) const
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }

   virtual void integrateRHSDirichletBoundary3D( const std::vector< Point3D >& coordsElement,
                                                 const std::vector< Point3D >& coordsFacet,
                                                 const Point3D&                oppositeVertex,
                                                 const Point3D&                outwardNormal,
                                                 const DGBasisInfo&            basis,
                                                 int                           degree,
                                                 MatrixXr&                     elMat ) const
   {
      WALBERLA_UNUSED( coordsElement );
      WALBERLA_UNUSED( coordsFacet );
      WALBERLA_UNUSED( oppositeVertex );
      WALBERLA_UNUSED( outwardNormal );
      WALBERLA_UNUSED( basis );
      WALBERLA_UNUSED( degree );
      WALBERLA_UNUSED( elMat );

      // Does nothing.
   }
};

} // namespace dg
} // namespace hyteg