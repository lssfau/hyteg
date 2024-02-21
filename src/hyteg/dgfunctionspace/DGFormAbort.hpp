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

using walberla::real_c;

/// \brief Helper form that aborts in every method.
///
/// Useful as a placeholder if a form is required but the functionality is not implemented yet.
class DGFormAbort : public DGForm
{
 protected:
   void integrateVolume2D( const std::vector< Point3D >&, const DGBasisInfo&,
                           const DGBasisInfo&,
                           int,
                           int, MatrixXr& )
       const override
   {
      WALBERLA_ABORT( "Invalid call to DGForm. This is a placeholder." );
   }

   virtual void integrateVolume3D( const std::vector< Point3D >&, const DGBasisInfo&,
                                   const DGBasisInfo&,
                                   int,
                                   int, MatrixXr& )
       const override
   {
      WALBERLA_ABORT( "Invalid call to DGForm. This is a placeholder." );
   }

   virtual void integrateFacetInner2D( const std::vector< Point3D >&,
                                       const std::vector< Point3D >&,
                                       const Point3D&,
                                       const Point3D&,
                                       const DGBasisInfo&,
                                       const DGBasisInfo&,
                                       int,
                                       int,
                                       MatrixXr& ) const
   {
      WALBERLA_ABORT( "Invalid call to DGForm. This is a placeholder." );
   }

   virtual void integrateFacetInner3D( const std::vector< Point3D >&,
                                       const std::vector< Point3D >&,
                                       const Point3D&,
                                       const Point3D&,
                                       const DGBasisInfo&,
                                       const DGBasisInfo&,
                                       int,
                                       int,
                                       MatrixXr& ) const
   {
      WALBERLA_ABORT( "Invalid call to DGForm. This is a placeholder." );
   }

   virtual void integrateFacetCoupling2D( const std::vector< Point3D >&,
                                          const std::vector< Point3D >&,
                                          const std::vector< Point3D >&,
                                          const Point3D&,
                                          const Point3D&,
                                          const Point3D&,
                                          const DGBasisInfo&,
                                          const DGBasisInfo&,
                                          int,
                                          int,
                                          MatrixXr& ) const
   {}

   virtual void integrateFacetCoupling3D( const std::vector< Point3D >&,
                                          const std::vector< Point3D >&,
                                          const std::vector< Point3D >&,
                                          const Point3D&,
                                          const Point3D&,
                                          const Point3D&,
                                          const DGBasisInfo&,
                                          const DGBasisInfo&,
                                          int,
                                          int,
                                          MatrixXr& ) const
   {
      WALBERLA_ABORT( "Invalid call to DGForm. This is a placeholder." );
   }

   virtual void integrateFacetDirichletBoundary2D( const std::vector< Point3D >&,
                                                   const std::vector< Point3D >&,
                                                   const Point3D&,
                                                   const Point3D&,
                                                   const DGBasisInfo&,
                                                   const DGBasisInfo&,
                                                   int,
                                                   int,
                                                   MatrixXr& ) const
   {
      WALBERLA_ABORT( "Invalid call to DGForm. This is a placeholder." );
   }

   virtual void integrateFacetDirichletBoundary3D( const std::vector< Point3D >&,
                                                   const std::vector< Point3D >&,
                                                   const Point3D&,
                                                   const Point3D&,
                                                   const DGBasisInfo&,
                                                   const DGBasisInfo&,
                                                   int,
                                                   int,
                                                   MatrixXr& ) const
   {
      WALBERLA_ABORT( "Invalid call to DGForm. This is a placeholder." );
   }

   virtual void integrateRHSDirichletBoundary2D( const std::vector< Point3D >&,
                                                 const std::vector< Point3D >&,
                                                 const Point3D&,
                                                 const Point3D&,
                                                 const DGBasisInfo&,
                                                 int,
                                                 MatrixXr& ) const
   {
      WALBERLA_ABORT( "Invalid call to DGForm. This is a placeholder." );
   }

   virtual void integrateRHSDirichletBoundary3D( const std::vector< Point3D >&,
                                                 const std::vector< Point3D >&,
                                                 const Point3D&,
                                                 const Point3D&,
                                                 const DGBasisInfo&,
                                                 int,
                                                 MatrixXr& ) const
   {
      WALBERLA_ABORT( "Invalid call to DGForm. This is a placeholder." );
   }
};

} // namespace dg
} // namespace hyteg