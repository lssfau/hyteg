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

class P1Form_div_1 : public P1FormHyTeG {
public:
  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
    Point3D x_hat({0.333333333333333, 0.333333333333333});
    Point3D x_tilde({0.333333333333333*coords[0][0] + 0.333333333333333*coords[1][0] + 0.333333333333333*coords[2][0], 0.333333333333333*coords[0][1] + 0.333333333333333*coords[1][1] + 0.333333333333333*coords[2][1]});
    Matrix2r DFinv;
    geometryMap_->evalDFinv(x_tilde, DFinv);
    real_t tmp0 = -coords[2][0];
    real_t tmp1 = -coords[2][1];
    real_t tmp2 = coords[0][0] - coords[1][0];
    real_t tmp3 = coords[0][1] + tmp1;
    real_t tmp4 = coords[0][0] + tmp0;
    real_t tmp5 = coords[0][1] - coords[1][1];
    real_t tmp6 = 1.0/(tmp2*tmp3 - tmp4*tmp5);
    real_t tmp7 = 0.166666666666667*tmp6/(fabs(tmp6)*fabs(DFinv(0,0)*DFinv(1,1) - DFinv(0,1)*DFinv(1,0)));
    out[0] = tmp7*(-DFinv(0,0)*(coords[1][1] + tmp1) + DFinv(1,0)*(coords[1][0] + tmp0));
    out[1] = tmp7*(DFinv(0,0)*tmp3 - DFinv(1,0)*tmp4);
    out[2] = tmp7*(-DFinv(0,0)*tmp5 + DFinv(1,0)*tmp2);
  }
};

class P1Form_div_2 : public P1FormHyTeG {
public:
  void integrate(const std::array<Point3D,3>& coords, Point3D& out) const
  {
    Point3D x_hat({0.333333333333333, 0.333333333333333});
    Point3D x_tilde({0.333333333333333*coords[0][0] + 0.333333333333333*coords[1][0] + 0.333333333333333*coords[2][0], 0.333333333333333*coords[0][1] + 0.333333333333333*coords[1][1] + 0.333333333333333*coords[2][1]});
    Matrix2r DFinv;
    geometryMap_->evalDFinv(x_tilde, DFinv);
    real_t tmp0 = -coords[2][0];
    real_t tmp1 = -coords[2][1];
    real_t tmp2 = coords[0][0] - coords[1][0];
    real_t tmp3 = coords[0][1] + tmp1;
    real_t tmp4 = coords[0][0] + tmp0;
    real_t tmp5 = coords[0][1] - coords[1][1];
    real_t tmp6 = 1.0/(tmp2*tmp3 - tmp4*tmp5);
    real_t tmp7 = 0.166666666666667*tmp6/(fabs(tmp6)*fabs(DFinv(0,0)*DFinv(1,1) - DFinv(0,1)*DFinv(1,0)));
    out[0] = tmp7*(-DFinv(0,1)*(coords[1][1] + tmp1) + DFinv(1,1)*(coords[1][0] + tmp0));
    out[1] = tmp7*(DFinv(0,1)*tmp3 - DFinv(1,1)*tmp4);
    out[2] = tmp7*(-DFinv(0,1)*tmp5 + DFinv(1,1)*tmp2);
  }
};

}
