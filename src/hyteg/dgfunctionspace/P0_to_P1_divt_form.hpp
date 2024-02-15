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
#include "hyteg/dgfunctionspace/DGForm2D.hpp"
#include "hyteg/dgfunctionspace/DGFormVolume.hpp"
#include "hyteg/types/Matrix.hpp"
#include "hyteg/types/PointND.hpp"

#include "Eigen/Eigen"

namespace hyteg {
namespace dg {

class p0_to_p1_divt_0_affine_q0 : public DGFormVolume
{
 protected:
   void integrateVolume2D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                     elMat ) const override;

   void integrateVolume3D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                     elMat ) const override;
};

class p0_to_p1_divt_1_affine_q0 : public DGFormVolume
{
 protected:
   void integrateVolume2D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                     elMat ) const override;

   void integrateVolume3D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                     elMat ) const override;
};

class p0_to_p1_divt_2_affine_q0 : public DGFormVolume
{
 protected:
   void integrateVolume2D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                     elMat ) const override;

   void integrateVolume3D( const std::vector< Point3D >& coords,
                           const DGBasisInfo&                                       trialBasis,
                           const DGBasisInfo&                                       testBasis,
                           int                                                      trialDegree,
                           int                                                      testDegree,
                           MatrixXr&                     elMat ) const override;
};

} // namespace dg
} // namespace hyteg