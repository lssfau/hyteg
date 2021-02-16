/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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

#include "hyteg/forms/form_hyteg_base/P1FormHyTeG.hpp"
#include "hyteg/geometry/GeometryMap.hpp"

namespace hyteg {

class P1Form_laplace : public P1FormHyTeG {
public:
  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
    Point3D x_hat({0.333333333333333, 0.333333333333333});
    Point3D x_tilde({0.333333333333333*coords[0][0] + 0.333333333333333*coords[1][0] + 0.333333333333333*coords[2][0], 0.333333333333333*coords[0][1] + 0.333333333333333*coords[1][1] + 0.333333333333333*coords[2][1]});
    Matrix2r DFinv;
    geometryMap_->evalDFinv(x_tilde, DFinv);
    real_t tmp0 = coords[0][1] - coords[1][1];
    real_t tmp1 = coords[0][0] - coords[1][0];
    real_t tmp2 = DFinv(0,0)*tmp0 - DFinv(1,0)*tmp1;
    real_t tmp3 = coords[0][0] - coords[2][0];
    real_t tmp4 = DFinv(1,0)*tmp3;
    real_t tmp5 = coords[0][1] - coords[2][1];
    real_t tmp6 = DFinv(0,0)*tmp5;
    real_t tmp7 = tmp2 + tmp4 - tmp6;
    real_t tmp8 = DFinv(0,1)*tmp0 - DFinv(1,1)*tmp1;
    real_t tmp9 = DFinv(1,1)*tmp3;
    real_t tmp10 = DFinv(0,1)*tmp5;
    real_t tmp11 = -tmp10 + tmp8 + tmp9;
    real_t tmp12 = 1.0/fabs(DFinv(0,0)*DFinv(1,1) - DFinv(0,1)*DFinv(1,0));
    real_t tmp13 = -tmp0*tmp3 + tmp1*tmp5;
    real_t tmp14 = pow(tmp13, -2);
    real_t tmp15 = 1.0/fabs(1.0/tmp13);
    real_t tmp16 = 0.5*tmp12*tmp14*tmp15;
    out[0] = tmp16*(pow(tmp11, 2) + pow(tmp7, 2));
    out[1] = tmp16*(tmp11*(tmp10 - tmp9) + tmp7*(-tmp4 + tmp6));
    out[2] = -tmp12*tmp14*tmp15*(0.5*tmp11*tmp8 + 0.5*tmp2*tmp7);
  }

  void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const
  {
    WALBERLA_ABORT( "P1Form_laplace not implemented for 3D" )
  }
};

}
