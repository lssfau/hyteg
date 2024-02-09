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

/// \brief Helper class to derive from if only 2D forms are implemented.
class DGForm2D : public DGForm
{
 protected:
   virtual void
       integrateVolume3D( const std::vector< Point3D >&, const DGBasisInfo&,
                                   const DGBasisInfo&,
                                   int,
                                   int, MatrixXr& ) const
   {
      WALBERLA_ABORT( "DGForm not implemented in 3D." );
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
      WALBERLA_ABORT( "DGForm not implemented in 3D." );
   }

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
      WALBERLA_ABORT( "DGForm not implemented in 3D." );
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
      WALBERLA_ABORT( "DGForm not implemented in 3D." );
   }

   virtual void integrateRHSDirichletBoundary3D( const std::vector< Point3D >&,
                                                 const std::vector< Point3D >&,
                                                 const Point3D&,
                                                 const Point3D&,
                                                 const DGBasisInfo&,
                                                 int,
                                                 MatrixXr& ) const
   {
      WALBERLA_ABORT( "DGForm not implemented in 3D." );
   }
};

} // namespace dg
} // namespace hyteg