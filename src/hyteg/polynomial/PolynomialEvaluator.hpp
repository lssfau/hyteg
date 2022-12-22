/*
 * Copyright (c) 2017-2021 Daniel Drzisga, Dominik Thoennes, Benjamin Mann
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

#include "Polynomial.hpp"

namespace hyteg {

namespace polynomialevaluator{
  template<uint_t Degree>
  inline real_t setStartX(real_t x, real_t h, const Polynomial1D<MonomialBasis1D>& poly1_, std::vector<real_t>& deltas) {
    static_assert(Degree <= 12, "Polynomial2DEvaluator not implemented for degree larger than 12");
    if (Degree == 0) {
      deltas[0] = poly1_.getCoefficient(0);
    }
    if (Degree == 1) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x;
      deltas[1] = h*poly1_.getCoefficient(1);
    }
    if (Degree == 2) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*x*x;
      deltas[1] = 2*h*poly1_.getCoefficient(2)*x + h*(h*poly1_.getCoefficient(2) + poly1_.getCoefficient(1));
      deltas[2] = 2*h*h*poly1_.getCoefficient(2);
    }
    if (Degree == 3) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x;
      deltas[1] = h*(h*(h*poly1_.getCoefficient(3) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(3*h*poly1_.getCoefficient(3)*x + h*(3*h*poly1_.getCoefficient(3) + 2*poly1_.getCoefficient(2)));
      deltas[2] = 6*h*h*poly1_.getCoefficient(3)*x + h*h*(6*h*poly1_.getCoefficient(3) + 2*poly1_.getCoefficient(2));
      deltas[3] = 6*h*h*h*poly1_.getCoefficient(3);
    }
    if (Degree == 4) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x + poly1_.getCoefficient(4)*x*x*x*x;
      deltas[1] = h*(h*(h*(h*poly1_.getCoefficient(4) + poly1_.getCoefficient(3)) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(h*(h*(4*h*poly1_.getCoefficient(4) + 3*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(4*h*poly1_.getCoefficient(4)*x + h*(6*h*poly1_.getCoefficient(4) + 3*poly1_.getCoefficient(3))));
      deltas[2] = h*h*(h*(14*h*poly1_.getCoefficient(4) + 6*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(12*h*h*poly1_.getCoefficient(4)*x + h*h*(24*h*poly1_.getCoefficient(4) + 6*poly1_.getCoefficient(3)));
      deltas[3] = 24*h*h*h*poly1_.getCoefficient(4)*x + h*h*h*(36*h*poly1_.getCoefficient(4) + 6*poly1_.getCoefficient(3));
      deltas[4] = 24*h*h*h*h*poly1_.getCoefficient(4);
    }
    if (Degree == 5) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x + poly1_.getCoefficient(4)*x*x*x*x + poly1_.getCoefficient(5)*x*x*x*x*x;
      deltas[1] = h*(h*(h*(h*(h*poly1_.getCoefficient(5) + poly1_.getCoefficient(4)) + poly1_.getCoefficient(3)) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(h*(h*(h*(5*h*poly1_.getCoefficient(5) + 4*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*(h*(10*h*poly1_.getCoefficient(5) + 6*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + x*(5*h*poly1_.getCoefficient(5)*x + h*(10*h*poly1_.getCoefficient(5) + 4*poly1_.getCoefficient(4)))));
      deltas[2] = h*h*(h*(h*(30*h*poly1_.getCoefficient(5) + 14*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*h*(h*(70*h*poly1_.getCoefficient(5) + 24*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(20*h*h*poly1_.getCoefficient(5)*x + h*h*(60*h*poly1_.getCoefficient(5) + 12*poly1_.getCoefficient(4))));
      deltas[3] = h*h*h*(h*(150*h*poly1_.getCoefficient(5) + 36*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(60*h*h*h*poly1_.getCoefficient(5)*x + h*h*h*(180*h*poly1_.getCoefficient(5) + 24*poly1_.getCoefficient(4)));
      deltas[4] = 120*h*h*h*h*poly1_.getCoefficient(5)*x + h*h*h*h*(240*h*poly1_.getCoefficient(5) + 24*poly1_.getCoefficient(4));
      deltas[5] = 120*h*h*h*h*h*poly1_.getCoefficient(5);
    }
    if (Degree == 6) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x + poly1_.getCoefficient(4)*x*x*x*x + poly1_.getCoefficient(5)*x*x*x*x*x + poly1_.getCoefficient(6)*x*x*x*x*x*x;
      deltas[1] = h*(h*(h*(h*(h*(h*poly1_.getCoefficient(6) + poly1_.getCoefficient(5)) + poly1_.getCoefficient(4)) + poly1_.getCoefficient(3)) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(h*(h*(h*(h*(6*h*poly1_.getCoefficient(6) + 5*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*(h*(h*(15*h*poly1_.getCoefficient(6) + 10*poly1_.getCoefficient(5)) + 6*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + x*(h*(h*(20*h*poly1_.getCoefficient(6) + 10*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + x*(6*h*poly1_.getCoefficient(6)*x + h*(15*h*poly1_.getCoefficient(6) + 5*poly1_.getCoefficient(5))))));
      deltas[2] = h*h*(h*(h*(h*(62*h*poly1_.getCoefficient(6) + 30*poly1_.getCoefficient(5)) + 14*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*h*(h*(h*(180*h*poly1_.getCoefficient(6) + 70*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*(h*(210*h*poly1_.getCoefficient(6) + 60*poly1_.getCoefficient(5)) + 12*poly1_.getCoefficient(4)) + x*(30*h*h*poly1_.getCoefficient(6)*x + h*h*(120*h*poly1_.getCoefficient(6) + 20*poly1_.getCoefficient(5)))));
      deltas[3] = h*h*h*(h*(h*(540*h*poly1_.getCoefficient(6) + 150*poly1_.getCoefficient(5)) + 36*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*h*(h*(900*h*poly1_.getCoefficient(6) + 180*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(120*h*h*h*poly1_.getCoefficient(6)*x + h*h*h*(540*h*poly1_.getCoefficient(6) + 60*poly1_.getCoefficient(5))));
      deltas[4] = h*h*h*h*(h*(1560*h*poly1_.getCoefficient(6) + 240*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(360*h*h*h*h*poly1_.getCoefficient(6)*x + h*h*h*h*(1440*h*poly1_.getCoefficient(6) + 120*poly1_.getCoefficient(5)));
      deltas[5] = 720*h*h*h*h*h*poly1_.getCoefficient(6)*x + h*h*h*h*h*(1800*h*poly1_.getCoefficient(6) + 120*poly1_.getCoefficient(5));
      deltas[6] = 720*h*h*h*h*h*h*poly1_.getCoefficient(6);
    }
    if (Degree == 7) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x + poly1_.getCoefficient(4)*x*x*x*x + poly1_.getCoefficient(5)*x*x*x*x*x + poly1_.getCoefficient(6)*x*x*x*x*x*x + poly1_.getCoefficient(7)*x*x*x*x*x*x*x;
      deltas[1] = h*(h*(h*(h*(h*(h*(h*poly1_.getCoefficient(7) + poly1_.getCoefficient(6)) + poly1_.getCoefficient(5)) + poly1_.getCoefficient(4)) + poly1_.getCoefficient(3)) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(h*(h*(h*(h*(h*(7*h*poly1_.getCoefficient(7) + 6*poly1_.getCoefficient(6)) + 5*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*(h*(h*(h*(21*h*poly1_.getCoefficient(7) + 15*poly1_.getCoefficient(6)) + 10*poly1_.getCoefficient(5)) + 6*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + x*(h*(h*(h*(35*h*poly1_.getCoefficient(7) + 20*poly1_.getCoefficient(6)) + 10*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + x*(h*(h*(35*h*poly1_.getCoefficient(7) + 15*poly1_.getCoefficient(6)) + 5*poly1_.getCoefficient(5)) + x*(7*h*poly1_.getCoefficient(7)*x + h*(21*h*poly1_.getCoefficient(7) + 6*poly1_.getCoefficient(6)))))));
      deltas[2] = h*h*(h*(h*(h*(h*(126*h*poly1_.getCoefficient(7) + 62*poly1_.getCoefficient(6)) + 30*poly1_.getCoefficient(5)) + 14*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*h*(h*(h*(h*(434*h*poly1_.getCoefficient(7) + 180*poly1_.getCoefficient(6)) + 70*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*(h*(h*(630*h*poly1_.getCoefficient(7) + 210*poly1_.getCoefficient(6)) + 60*poly1_.getCoefficient(5)) + 12*poly1_.getCoefficient(4)) + x*(h*h*(h*(490*h*poly1_.getCoefficient(7) + 120*poly1_.getCoefficient(6)) + 20*poly1_.getCoefficient(5)) + x*(42*h*h*poly1_.getCoefficient(7)*x + h*h*(210*h*poly1_.getCoefficient(7) + 30*poly1_.getCoefficient(6))))));
      deltas[3] = h*h*h*(h*(h*(h*(1806*h*poly1_.getCoefficient(7) + 540*poly1_.getCoefficient(6)) + 150*poly1_.getCoefficient(5)) + 36*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*h*(h*(h*(3780*h*poly1_.getCoefficient(7) + 900*poly1_.getCoefficient(6)) + 180*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(h*h*h*(h*(3150*h*poly1_.getCoefficient(7) + 540*poly1_.getCoefficient(6)) + 60*poly1_.getCoefficient(5)) + x*(210*h*h*h*poly1_.getCoefficient(7)*x + h*h*h*(1260*h*poly1_.getCoefficient(7) + 120*poly1_.getCoefficient(6)))));
      deltas[4] = h*h*h*h*(h*(h*(8400*h*poly1_.getCoefficient(7) + 1560*poly1_.getCoefficient(6)) + 240*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(h*h*h*h*(h*(10920*h*poly1_.getCoefficient(7) + 1440*poly1_.getCoefficient(6)) + 120*poly1_.getCoefficient(5)) + x*(840*h*h*h*h*poly1_.getCoefficient(7)*x + h*h*h*h*(5040*h*poly1_.getCoefficient(7) + 360*poly1_.getCoefficient(6))));
      deltas[5] = h*h*h*h*h*(h*(16800*h*poly1_.getCoefficient(7) + 1800*poly1_.getCoefficient(6)) + 120*poly1_.getCoefficient(5)) + x*(2520*h*h*h*h*h*poly1_.getCoefficient(7)*x + h*h*h*h*h*(12600*h*poly1_.getCoefficient(7) + 720*poly1_.getCoefficient(6)));
      deltas[6] = 5040*h*h*h*h*h*h*poly1_.getCoefficient(7)*x + h*h*h*h*h*h*(15120*h*poly1_.getCoefficient(7) + 720*poly1_.getCoefficient(6));
      deltas[7] = 5040*h*h*h*h*h*h*h*poly1_.getCoefficient(7);
    }
    if (Degree == 8) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x + poly1_.getCoefficient(4)*x*x*x*x + poly1_.getCoefficient(5)*x*x*x*x*x + poly1_.getCoefficient(6)*x*x*x*x*x*x + poly1_.getCoefficient(7)*x*x*x*x*x*x*x + poly1_.getCoefficient(8)*x*x*x*x*x*x*x*x;
      deltas[1] = h*(h*(h*(h*(h*(h*(h*(h*poly1_.getCoefficient(8) + poly1_.getCoefficient(7)) + poly1_.getCoefficient(6)) + poly1_.getCoefficient(5)) + poly1_.getCoefficient(4)) + poly1_.getCoefficient(3)) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(h*(h*(h*(h*(h*(h*(8*h*poly1_.getCoefficient(8) + 7*poly1_.getCoefficient(7)) + 6*poly1_.getCoefficient(6)) + 5*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*(h*(h*(h*(h*(28*h*poly1_.getCoefficient(8) + 21*poly1_.getCoefficient(7)) + 15*poly1_.getCoefficient(6)) + 10*poly1_.getCoefficient(5)) + 6*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + x*(h*(h*(h*(h*(56*h*poly1_.getCoefficient(8) + 35*poly1_.getCoefficient(7)) + 20*poly1_.getCoefficient(6)) + 10*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + x*(h*(h*(h*(70*h*poly1_.getCoefficient(8) + 35*poly1_.getCoefficient(7)) + 15*poly1_.getCoefficient(6)) + 5*poly1_.getCoefficient(5)) + x*(h*(h*(56*h*poly1_.getCoefficient(8) + 21*poly1_.getCoefficient(7)) + 6*poly1_.getCoefficient(6)) + x*(8*h*poly1_.getCoefficient(8)*x + h*(28*h*poly1_.getCoefficient(8) + 7*poly1_.getCoefficient(7))))))));
      deltas[2] = h*h*(h*(h*(h*(h*(h*(254*h*poly1_.getCoefficient(8) + 126*poly1_.getCoefficient(7)) + 62*poly1_.getCoefficient(6)) + 30*poly1_.getCoefficient(5)) + 14*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*h*(h*(h*(h*(h*(1008*h*poly1_.getCoefficient(8) + 434*poly1_.getCoefficient(7)) + 180*poly1_.getCoefficient(6)) + 70*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*(h*(h*(h*(1736*h*poly1_.getCoefficient(8) + 630*poly1_.getCoefficient(7)) + 210*poly1_.getCoefficient(6)) + 60*poly1_.getCoefficient(5)) + 12*poly1_.getCoefficient(4)) + x*(h*h*(h*(h*(1680*h*poly1_.getCoefficient(8) + 490*poly1_.getCoefficient(7)) + 120*poly1_.getCoefficient(6)) + 20*poly1_.getCoefficient(5)) + x*(h*h*(h*(980*h*poly1_.getCoefficient(8) + 210*poly1_.getCoefficient(7)) + 30*poly1_.getCoefficient(6)) + x*(56*h*h*poly1_.getCoefficient(8)*x + h*h*(336*h*poly1_.getCoefficient(8) + 42*poly1_.getCoefficient(7)))))));
      deltas[3] = h*h*h*(h*(h*(h*(h*(5796*h*poly1_.getCoefficient(8) + 1806*poly1_.getCoefficient(7)) + 540*poly1_.getCoefficient(6)) + 150*poly1_.getCoefficient(5)) + 36*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*h*(h*(h*(h*(14448*h*poly1_.getCoefficient(8) + 3780*poly1_.getCoefficient(7)) + 900*poly1_.getCoefficient(6)) + 180*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(h*h*h*(h*(h*(15120*h*poly1_.getCoefficient(8) + 3150*poly1_.getCoefficient(7)) + 540*poly1_.getCoefficient(6)) + 60*poly1_.getCoefficient(5)) + x*(h*h*h*(h*(8400*h*poly1_.getCoefficient(8) + 1260*poly1_.getCoefficient(7)) + 120*poly1_.getCoefficient(6)) + x*(336*h*h*h*poly1_.getCoefficient(8)*x + h*h*h*(2520*h*poly1_.getCoefficient(8) + 210*poly1_.getCoefficient(7))))));
      deltas[4] = h*h*h*h*(h*(h*(h*(40824*h*poly1_.getCoefficient(8) + 8400*poly1_.getCoefficient(7)) + 1560*poly1_.getCoefficient(6)) + 240*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(h*h*h*h*(h*(h*(67200*h*poly1_.getCoefficient(8) + 10920*poly1_.getCoefficient(7)) + 1440*poly1_.getCoefficient(6)) + 120*poly1_.getCoefficient(5)) + x*(h*h*h*h*(h*(43680*h*poly1_.getCoefficient(8) + 5040*poly1_.getCoefficient(7)) + 360*poly1_.getCoefficient(6)) + x*(1680*h*h*h*h*poly1_.getCoefficient(8)*x + h*h*h*h*(13440*h*poly1_.getCoefficient(8) + 840*poly1_.getCoefficient(7)))));
      deltas[5] = h*h*h*h*h*(h*(h*(126000*h*poly1_.getCoefficient(8) + 16800*poly1_.getCoefficient(7)) + 1800*poly1_.getCoefficient(6)) + 120*poly1_.getCoefficient(5)) + x*(h*h*h*h*h*(h*(134400*h*poly1_.getCoefficient(8) + 12600*poly1_.getCoefficient(7)) + 720*poly1_.getCoefficient(6)) + x*(6720*h*h*h*h*h*poly1_.getCoefficient(8)*x + h*h*h*h*h*(50400*h*poly1_.getCoefficient(8) + 2520*poly1_.getCoefficient(7))));
      deltas[6] = h*h*h*h*h*h*(h*(191520*h*poly1_.getCoefficient(8) + 15120*poly1_.getCoefficient(7)) + 720*poly1_.getCoefficient(6)) + x*(20160*h*h*h*h*h*h*poly1_.getCoefficient(8)*x + h*h*h*h*h*h*(120960*h*poly1_.getCoefficient(8) + 5040*poly1_.getCoefficient(7)));
      deltas[7] = 40320*h*h*h*h*h*h*h*poly1_.getCoefficient(8)*x + h*h*h*h*h*h*h*(141120*h*poly1_.getCoefficient(8) + 5040*poly1_.getCoefficient(7));
      deltas[8] = 40320*h*h*h*h*h*h*h*h*poly1_.getCoefficient(8);
    }
    if (Degree == 9) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x + poly1_.getCoefficient(4)*x*x*x*x + poly1_.getCoefficient(5)*x*x*x*x*x + poly1_.getCoefficient(6)*x*x*x*x*x*x + poly1_.getCoefficient(7)*x*x*x*x*x*x*x + poly1_.getCoefficient(8)*x*x*x*x*x*x*x*x + poly1_.getCoefficient(9)*x*x*x*x*x*x*x*x*x;
      deltas[1] = h*(h*(h*(h*(h*(h*(h*(h*(h*poly1_.getCoefficient(9) + poly1_.getCoefficient(8)) + poly1_.getCoefficient(7)) + poly1_.getCoefficient(6)) + poly1_.getCoefficient(5)) + poly1_.getCoefficient(4)) + poly1_.getCoefficient(3)) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(h*(h*(h*(h*(h*(h*(h*(9*h*poly1_.getCoefficient(9) + 8*poly1_.getCoefficient(8)) + 7*poly1_.getCoefficient(7)) + 6*poly1_.getCoefficient(6)) + 5*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*(h*(h*(h*(h*(h*(36*h*poly1_.getCoefficient(9) + 28*poly1_.getCoefficient(8)) + 21*poly1_.getCoefficient(7)) + 15*poly1_.getCoefficient(6)) + 10*poly1_.getCoefficient(5)) + 6*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + x*(h*(h*(h*(h*(h*(84*h*poly1_.getCoefficient(9) + 56*poly1_.getCoefficient(8)) + 35*poly1_.getCoefficient(7)) + 20*poly1_.getCoefficient(6)) + 10*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + x*(h*(h*(h*(h*(126*h*poly1_.getCoefficient(9) + 70*poly1_.getCoefficient(8)) + 35*poly1_.getCoefficient(7)) + 15*poly1_.getCoefficient(6)) + 5*poly1_.getCoefficient(5)) + x*(h*(h*(h*(126*h*poly1_.getCoefficient(9) + 56*poly1_.getCoefficient(8)) + 21*poly1_.getCoefficient(7)) + 6*poly1_.getCoefficient(6)) + x*(h*(h*(84*h*poly1_.getCoefficient(9) + 28*poly1_.getCoefficient(8)) + 7*poly1_.getCoefficient(7)) + x*(9*h*poly1_.getCoefficient(9)*x + h*(36*h*poly1_.getCoefficient(9) + 8*poly1_.getCoefficient(8)))))))));
      deltas[2] = h*h*(h*(h*(h*(h*(h*(h*(510*h*poly1_.getCoefficient(9) + 254*poly1_.getCoefficient(8)) + 126*poly1_.getCoefficient(7)) + 62*poly1_.getCoefficient(6)) + 30*poly1_.getCoefficient(5)) + 14*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*h*(h*(h*(h*(h*(h*(2286*h*poly1_.getCoefficient(9) + 1008*poly1_.getCoefficient(8)) + 434*poly1_.getCoefficient(7)) + 180*poly1_.getCoefficient(6)) + 70*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*(h*(h*(h*(h*(4536*h*poly1_.getCoefficient(9) + 1736*poly1_.getCoefficient(8)) + 630*poly1_.getCoefficient(7)) + 210*poly1_.getCoefficient(6)) + 60*poly1_.getCoefficient(5)) + 12*poly1_.getCoefficient(4)) + x*(h*h*(h*(h*(h*(5208*h*poly1_.getCoefficient(9) + 1680*poly1_.getCoefficient(8)) + 490*poly1_.getCoefficient(7)) + 120*poly1_.getCoefficient(6)) + 20*poly1_.getCoefficient(5)) + x*(h*h*(h*(h*(3780*h*poly1_.getCoefficient(9) + 980*poly1_.getCoefficient(8)) + 210*poly1_.getCoefficient(7)) + 30*poly1_.getCoefficient(6)) + x*(h*h*(h*(1764*h*poly1_.getCoefficient(9) + 336*poly1_.getCoefficient(8)) + 42*poly1_.getCoefficient(7)) + x*(72*h*h*poly1_.getCoefficient(9)*x + h*h*(504*h*poly1_.getCoefficient(9) + 56*poly1_.getCoefficient(8))))))));
      deltas[3] = h*h*h*(h*(h*(h*(h*(h*(18150*h*poly1_.getCoefficient(9) + 5796*poly1_.getCoefficient(8)) + 1806*poly1_.getCoefficient(7)) + 540*poly1_.getCoefficient(6)) + 150*poly1_.getCoefficient(5)) + 36*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*h*(h*(h*(h*(h*(52164*h*poly1_.getCoefficient(9) + 14448*poly1_.getCoefficient(8)) + 3780*poly1_.getCoefficient(7)) + 900*poly1_.getCoefficient(6)) + 180*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(h*h*h*(h*(h*(h*(65016*h*poly1_.getCoefficient(9) + 15120*poly1_.getCoefficient(8)) + 3150*poly1_.getCoefficient(7)) + 540*poly1_.getCoefficient(6)) + 60*poly1_.getCoefficient(5)) + x*(h*h*h*(h*(h*(45360*h*poly1_.getCoefficient(9) + 8400*poly1_.getCoefficient(8)) + 1260*poly1_.getCoefficient(7)) + 120*poly1_.getCoefficient(6)) + x*(h*h*h*(h*(18900*h*poly1_.getCoefficient(9) + 2520*poly1_.getCoefficient(8)) + 210*poly1_.getCoefficient(7)) + x*(504*h*h*h*poly1_.getCoefficient(9)*x + h*h*h*(4536*h*poly1_.getCoefficient(9) + 336*poly1_.getCoefficient(8)))))));
      deltas[4] = h*h*h*h*(h*(h*(h*(h*(186480*h*poly1_.getCoefficient(9) + 40824*poly1_.getCoefficient(8)) + 8400*poly1_.getCoefficient(7)) + 1560*poly1_.getCoefficient(6)) + 240*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(h*h*h*h*(h*(h*(h*(367416*h*poly1_.getCoefficient(9) + 67200*poly1_.getCoefficient(8)) + 10920*poly1_.getCoefficient(7)) + 1440*poly1_.getCoefficient(6)) + 120*poly1_.getCoefficient(5)) + x*(h*h*h*h*(h*(h*(302400*h*poly1_.getCoefficient(9) + 43680*poly1_.getCoefficient(8)) + 5040*poly1_.getCoefficient(7)) + 360*poly1_.getCoefficient(6)) + x*(h*h*h*h*(h*(131040*h*poly1_.getCoefficient(9) + 13440*poly1_.getCoefficient(8)) + 840*poly1_.getCoefficient(7)) + x*(3024*h*h*h*h*poly1_.getCoefficient(9)*x + h*h*h*h*(30240*h*poly1_.getCoefficient(9) + 1680*poly1_.getCoefficient(8))))));
      deltas[5] = h*h*h*h*h*(h*(h*(h*(834120*h*poly1_.getCoefficient(9) + 126000*poly1_.getCoefficient(8)) + 16800*poly1_.getCoefficient(7)) + 1800*poly1_.getCoefficient(6)) + 120*poly1_.getCoefficient(5)) + x*(h*h*h*h*h*(h*(h*(1134000*h*poly1_.getCoefficient(9) + 134400*poly1_.getCoefficient(8)) + 12600*poly1_.getCoefficient(7)) + 720*poly1_.getCoefficient(6)) + x*(h*h*h*h*h*(h*(604800*h*poly1_.getCoefficient(9) + 50400*poly1_.getCoefficient(8)) + 2520*poly1_.getCoefficient(7)) + x*(15120*h*h*h*h*h*poly1_.getCoefficient(9)*x + h*h*h*h*h*(151200*h*poly1_.getCoefficient(9) + 6720*poly1_.getCoefficient(8)))));
      deltas[6] = h*h*h*h*h*h*(h*(h*(1905120*h*poly1_.getCoefficient(9) + 191520*poly1_.getCoefficient(8)) + 15120*poly1_.getCoefficient(7)) + 720*poly1_.getCoefficient(6)) + x*(h*h*h*h*h*h*(h*(1723680*h*poly1_.getCoefficient(9) + 120960*poly1_.getCoefficient(8)) + 5040*poly1_.getCoefficient(7)) + x*(60480*h*h*h*h*h*h*poly1_.getCoefficient(9)*x + h*h*h*h*h*h*(544320*h*poly1_.getCoefficient(9) + 20160*poly1_.getCoefficient(8))));
      deltas[7] = h*h*h*h*h*h*h*(h*(2328480*h*poly1_.getCoefficient(9) + 141120*poly1_.getCoefficient(8)) + 5040*poly1_.getCoefficient(7)) + x*(181440*h*h*h*h*h*h*h*poly1_.getCoefficient(9)*x + h*h*h*h*h*h*h*(1270080*h*poly1_.getCoefficient(9) + 40320*poly1_.getCoefficient(8)));
      deltas[8] = 362880*h*h*h*h*h*h*h*h*poly1_.getCoefficient(9)*x + h*h*h*h*h*h*h*h*(1451520*h*poly1_.getCoefficient(9) + 40320*poly1_.getCoefficient(8));
      deltas[9] = 362880*h*h*h*h*h*h*h*h*h*poly1_.getCoefficient(9);
    }
    if (Degree == 10) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(10)*x*x*x*x*x*x*x*x*x*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x + poly1_.getCoefficient(4)*x*x*x*x + poly1_.getCoefficient(5)*x*x*x*x*x + poly1_.getCoefficient(6)*x*x*x*x*x*x + poly1_.getCoefficient(7)*x*x*x*x*x*x*x + poly1_.getCoefficient(8)*x*x*x*x*x*x*x*x + poly1_.getCoefficient(9)*x*x*x*x*x*x*x*x*x;
      deltas[1] = h*(h*(h*(h*(h*(h*(h*(h*(h*(h*poly1_.getCoefficient(10) + poly1_.getCoefficient(9)) + poly1_.getCoefficient(8)) + poly1_.getCoefficient(7)) + poly1_.getCoefficient(6)) + poly1_.getCoefficient(5)) + poly1_.getCoefficient(4)) + poly1_.getCoefficient(3)) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(h*(h*(h*(h*(h*(h*(h*(h*(10*h*poly1_.getCoefficient(10) + 9*poly1_.getCoefficient(9)) + 8*poly1_.getCoefficient(8)) + 7*poly1_.getCoefficient(7)) + 6*poly1_.getCoefficient(6)) + 5*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*(h*(h*(h*(h*(h*(h*(45*h*poly1_.getCoefficient(10) + 36*poly1_.getCoefficient(9)) + 28*poly1_.getCoefficient(8)) + 21*poly1_.getCoefficient(7)) + 15*poly1_.getCoefficient(6)) + 10*poly1_.getCoefficient(5)) + 6*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + x*(h*(h*(h*(h*(h*(h*(120*h*poly1_.getCoefficient(10) + 84*poly1_.getCoefficient(9)) + 56*poly1_.getCoefficient(8)) + 35*poly1_.getCoefficient(7)) + 20*poly1_.getCoefficient(6)) + 10*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + x*(h*(h*(h*(h*(h*(210*h*poly1_.getCoefficient(10) + 126*poly1_.getCoefficient(9)) + 70*poly1_.getCoefficient(8)) + 35*poly1_.getCoefficient(7)) + 15*poly1_.getCoefficient(6)) + 5*poly1_.getCoefficient(5)) + x*(h*(h*(h*(h*(252*h*poly1_.getCoefficient(10) + 126*poly1_.getCoefficient(9)) + 56*poly1_.getCoefficient(8)) + 21*poly1_.getCoefficient(7)) + 6*poly1_.getCoefficient(6)) + x*(h*(h*(h*(210*h*poly1_.getCoefficient(10) + 84*poly1_.getCoefficient(9)) + 28*poly1_.getCoefficient(8)) + 7*poly1_.getCoefficient(7)) + x*(h*(h*(120*h*poly1_.getCoefficient(10) + 36*poly1_.getCoefficient(9)) + 8*poly1_.getCoefficient(8)) + x*(10*h*poly1_.getCoefficient(10)*x + h*(45*h*poly1_.getCoefficient(10) + 9*poly1_.getCoefficient(9))))))))));
      deltas[2] = h*h*(h*(h*(h*(h*(h*(h*(h*(1022*h*poly1_.getCoefficient(10) + 510*poly1_.getCoefficient(9)) + 254*poly1_.getCoefficient(8)) + 126*poly1_.getCoefficient(7)) + 62*poly1_.getCoefficient(6)) + 30*poly1_.getCoefficient(5)) + 14*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*h*(h*(h*(h*(h*(h*(h*(5100*h*poly1_.getCoefficient(10) + 2286*poly1_.getCoefficient(9)) + 1008*poly1_.getCoefficient(8)) + 434*poly1_.getCoefficient(7)) + 180*poly1_.getCoefficient(6)) + 70*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*(h*(h*(h*(h*(h*(11430*h*poly1_.getCoefficient(10) + 4536*poly1_.getCoefficient(9)) + 1736*poly1_.getCoefficient(8)) + 630*poly1_.getCoefficient(7)) + 210*poly1_.getCoefficient(6)) + 60*poly1_.getCoefficient(5)) + 12*poly1_.getCoefficient(4)) + x*(h*h*(h*(h*(h*(h*(15120*h*poly1_.getCoefficient(10) + 5208*poly1_.getCoefficient(9)) + 1680*poly1_.getCoefficient(8)) + 490*poly1_.getCoefficient(7)) + 120*poly1_.getCoefficient(6)) + 20*poly1_.getCoefficient(5)) + x*(h*h*(h*(h*(h*(13020*h*poly1_.getCoefficient(10) + 3780*poly1_.getCoefficient(9)) + 980*poly1_.getCoefficient(8)) + 210*poly1_.getCoefficient(7)) + 30*poly1_.getCoefficient(6)) + x*(h*h*(h*(h*(7560*h*poly1_.getCoefficient(10) + 1764*poly1_.getCoefficient(9)) + 336*poly1_.getCoefficient(8)) + 42*poly1_.getCoefficient(7)) + x*(h*h*(h*(2940*h*poly1_.getCoefficient(10) + 504*poly1_.getCoefficient(9)) + 56*poly1_.getCoefficient(8)) + x*(90*h*h*poly1_.getCoefficient(10)*x + h*h*(720*h*poly1_.getCoefficient(10) + 72*poly1_.getCoefficient(9)))))))));
      deltas[3] = h*h*h*(h*(h*(h*(h*(h*(h*(55980*h*poly1_.getCoefficient(10) + 18150*poly1_.getCoefficient(9)) + 5796*poly1_.getCoefficient(8)) + 1806*poly1_.getCoefficient(7)) + 540*poly1_.getCoefficient(6)) + 150*poly1_.getCoefficient(5)) + 36*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*h*(h*(h*(h*(h*(h*(181500*h*poly1_.getCoefficient(10) + 52164*poly1_.getCoefficient(9)) + 14448*poly1_.getCoefficient(8)) + 3780*poly1_.getCoefficient(7)) + 900*poly1_.getCoefficient(6)) + 180*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(h*h*h*(h*(h*(h*(h*(260820*h*poly1_.getCoefficient(10) + 65016*poly1_.getCoefficient(9)) + 15120*poly1_.getCoefficient(8)) + 3150*poly1_.getCoefficient(7)) + 540*poly1_.getCoefficient(6)) + 60*poly1_.getCoefficient(5)) + x*(h*h*h*(h*(h*(h*(216720*h*poly1_.getCoefficient(10) + 45360*poly1_.getCoefficient(9)) + 8400*poly1_.getCoefficient(8)) + 1260*poly1_.getCoefficient(7)) + 120*poly1_.getCoefficient(6)) + x*(h*h*h*(h*(h*(113400*h*poly1_.getCoefficient(10) + 18900*poly1_.getCoefficient(9)) + 2520*poly1_.getCoefficient(8)) + 210*poly1_.getCoefficient(7)) + x*(h*h*h*(h*(37800*h*poly1_.getCoefficient(10) + 4536*poly1_.getCoefficient(9)) + 336*poly1_.getCoefficient(8)) + x*(720*h*h*h*poly1_.getCoefficient(10)*x + h*h*h*(7560*h*poly1_.getCoefficient(10) + 504*poly1_.getCoefficient(9))))))));
      deltas[4] = h*h*h*h*(h*(h*(h*(h*(h*(818520*h*poly1_.getCoefficient(10) + 186480*poly1_.getCoefficient(9)) + 40824*poly1_.getCoefficient(8)) + 8400*poly1_.getCoefficient(7)) + 1560*poly1_.getCoefficient(6)) + 240*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(h*h*h*h*(h*(h*(h*(h*(1864800*h*poly1_.getCoefficient(10) + 367416*poly1_.getCoefficient(9)) + 67200*poly1_.getCoefficient(8)) + 10920*poly1_.getCoefficient(7)) + 1440*poly1_.getCoefficient(6)) + 120*poly1_.getCoefficient(5)) + x*(h*h*h*h*(h*(h*(h*(1837080*h*poly1_.getCoefficient(10) + 302400*poly1_.getCoefficient(9)) + 43680*poly1_.getCoefficient(8)) + 5040*poly1_.getCoefficient(7)) + 360*poly1_.getCoefficient(6)) + x*(h*h*h*h*(h*(h*(1008000*h*poly1_.getCoefficient(10) + 131040*poly1_.getCoefficient(9)) + 13440*poly1_.getCoefficient(8)) + 840*poly1_.getCoefficient(7)) + x*(h*h*h*h*(h*(327600*h*poly1_.getCoefficient(10) + 30240*poly1_.getCoefficient(9)) + 1680*poly1_.getCoefficient(8)) + x*(5040*h*h*h*h*poly1_.getCoefficient(10)*x + h*h*h*h*(60480*h*poly1_.getCoefficient(10) + 3024*poly1_.getCoefficient(9)))))));
      deltas[5] = h*h*h*h*h*(h*(h*(h*(h*(5103000*h*poly1_.getCoefficient(10) + 834120*poly1_.getCoefficient(9)) + 126000*poly1_.getCoefficient(8)) + 16800*poly1_.getCoefficient(7)) + 1800*poly1_.getCoefficient(6)) + 120*poly1_.getCoefficient(5)) + x*(h*h*h*h*h*(h*(h*(h*(8341200*h*poly1_.getCoefficient(10) + 1134000*poly1_.getCoefficient(9)) + 134400*poly1_.getCoefficient(8)) + 12600*poly1_.getCoefficient(7)) + 720*poly1_.getCoefficient(6)) + x*(h*h*h*h*h*(h*(h*(5670000*h*poly1_.getCoefficient(10) + 604800*poly1_.getCoefficient(9)) + 50400*poly1_.getCoefficient(8)) + 2520*poly1_.getCoefficient(7)) + x*(h*h*h*h*h*(h*(2016000*h*poly1_.getCoefficient(10) + 151200*poly1_.getCoefficient(9)) + 6720*poly1_.getCoefficient(8)) + x*(30240*h*h*h*h*h*poly1_.getCoefficient(10)*x + h*h*h*h*h*(378000*h*poly1_.getCoefficient(10) + 15120*poly1_.getCoefficient(9))))));
      deltas[6] = h*h*h*h*h*h*(h*(h*(h*(16435440*h*poly1_.getCoefficient(10) + 1905120*poly1_.getCoefficient(9)) + 191520*poly1_.getCoefficient(8)) + 15120*poly1_.getCoefficient(7)) + 720*poly1_.getCoefficient(6)) + x*(h*h*h*h*h*h*(h*(h*(19051200*h*poly1_.getCoefficient(10) + 1723680*poly1_.getCoefficient(9)) + 120960*poly1_.getCoefficient(8)) + 5040*poly1_.getCoefficient(7)) + x*(h*h*h*h*h*h*(h*(8618400*h*poly1_.getCoefficient(10) + 544320*poly1_.getCoefficient(9)) + 20160*poly1_.getCoefficient(8)) + x*(151200*h*h*h*h*h*h*poly1_.getCoefficient(10)*x + h*h*h*h*h*h*(1814400*h*poly1_.getCoefficient(10) + 60480*poly1_.getCoefficient(9)))));
      deltas[7] = h*h*h*h*h*h*h*(h*(h*(29635200*h*poly1_.getCoefficient(10) + 2328480*poly1_.getCoefficient(9)) + 141120*poly1_.getCoefficient(8)) + 5040*poly1_.getCoefficient(7)) + x*(h*h*h*h*h*h*h*(h*(23284800*h*poly1_.getCoefficient(10) + 1270080*poly1_.getCoefficient(9)) + 40320*poly1_.getCoefficient(8)) + x*(604800*h*h*h*h*h*h*h*poly1_.getCoefficient(10)*x + h*h*h*h*h*h*h*(6350400*h*poly1_.getCoefficient(10) + 181440*poly1_.getCoefficient(9))));
      deltas[8] = h*h*h*h*h*h*h*h*(h*(30240000*h*poly1_.getCoefficient(10) + 1451520*poly1_.getCoefficient(9)) + 40320*poly1_.getCoefficient(8)) + x*(1814400*h*h*h*h*h*h*h*h*poly1_.getCoefficient(10)*x + h*h*h*h*h*h*h*h*(14515200*h*poly1_.getCoefficient(10) + 362880*poly1_.getCoefficient(9)));
      deltas[9] = 3628800*h*h*h*h*h*h*h*h*h*poly1_.getCoefficient(10)*x + h*h*h*h*h*h*h*h*h*(16329600*h*poly1_.getCoefficient(10) + 362880*poly1_.getCoefficient(9));
      deltas[10] = 3628800*h*h*h*h*h*h*h*h*h*h*poly1_.getCoefficient(10);
    }
    if (Degree == 11) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(10)*x*x*x*x*x*x*x*x*x*x + poly1_.getCoefficient(11)*x*x*x*x*x*x*x*x*x*x*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x + poly1_.getCoefficient(4)*x*x*x*x + poly1_.getCoefficient(5)*x*x*x*x*x + poly1_.getCoefficient(6)*x*x*x*x*x*x + poly1_.getCoefficient(7)*x*x*x*x*x*x*x + poly1_.getCoefficient(8)*x*x*x*x*x*x*x*x + poly1_.getCoefficient(9)*x*x*x*x*x*x*x*x*x;
      deltas[1] = h*(h*(h*(h*(h*(h*(h*(h*(h*(h*(h*poly1_.getCoefficient(11) + poly1_.getCoefficient(10)) + poly1_.getCoefficient(9)) + poly1_.getCoefficient(8)) + poly1_.getCoefficient(7)) + poly1_.getCoefficient(6)) + poly1_.getCoefficient(5)) + poly1_.getCoefficient(4)) + poly1_.getCoefficient(3)) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(h*(h*(h*(h*(h*(h*(h*(h*(h*(11*h*poly1_.getCoefficient(11) + 10*poly1_.getCoefficient(10)) + 9*poly1_.getCoefficient(9)) + 8*poly1_.getCoefficient(8)) + 7*poly1_.getCoefficient(7)) + 6*poly1_.getCoefficient(6)) + 5*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*(h*(h*(h*(h*(h*(h*(h*(55*h*poly1_.getCoefficient(11) + 45*poly1_.getCoefficient(10)) + 36*poly1_.getCoefficient(9)) + 28*poly1_.getCoefficient(8)) + 21*poly1_.getCoefficient(7)) + 15*poly1_.getCoefficient(6)) + 10*poly1_.getCoefficient(5)) + 6*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + x*(h*(h*(h*(h*(h*(h*(h*(165*h*poly1_.getCoefficient(11) + 120*poly1_.getCoefficient(10)) + 84*poly1_.getCoefficient(9)) + 56*poly1_.getCoefficient(8)) + 35*poly1_.getCoefficient(7)) + 20*poly1_.getCoefficient(6)) + 10*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + x*(h*(h*(h*(h*(h*(h*(330*h*poly1_.getCoefficient(11) + 210*poly1_.getCoefficient(10)) + 126*poly1_.getCoefficient(9)) + 70*poly1_.getCoefficient(8)) + 35*poly1_.getCoefficient(7)) + 15*poly1_.getCoefficient(6)) + 5*poly1_.getCoefficient(5)) + x*(h*(h*(h*(h*(h*(462*h*poly1_.getCoefficient(11) + 252*poly1_.getCoefficient(10)) + 126*poly1_.getCoefficient(9)) + 56*poly1_.getCoefficient(8)) + 21*poly1_.getCoefficient(7)) + 6*poly1_.getCoefficient(6)) + x*(h*(h*(h*(h*(462*h*poly1_.getCoefficient(11) + 210*poly1_.getCoefficient(10)) + 84*poly1_.getCoefficient(9)) + 28*poly1_.getCoefficient(8)) + 7*poly1_.getCoefficient(7)) + x*(h*(h*(h*(330*h*poly1_.getCoefficient(11) + 120*poly1_.getCoefficient(10)) + 36*poly1_.getCoefficient(9)) + 8*poly1_.getCoefficient(8)) + x*(h*(h*(165*h*poly1_.getCoefficient(11) + 45*poly1_.getCoefficient(10)) + 9*poly1_.getCoefficient(9)) + x*(11*h*poly1_.getCoefficient(11)*x + h*(55*h*poly1_.getCoefficient(11) + 10*poly1_.getCoefficient(10)))))))))));
      deltas[2] = h*h*(h*(h*(h*(h*(h*(h*(h*(h*(2046*h*poly1_.getCoefficient(11) + 1022*poly1_.getCoefficient(10)) + 510*poly1_.getCoefficient(9)) + 254*poly1_.getCoefficient(8)) + 126*poly1_.getCoefficient(7)) + 62*poly1_.getCoefficient(6)) + 30*poly1_.getCoefficient(5)) + 14*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*h*(h*(h*(h*(h*(h*(h*(h*(11242*h*poly1_.getCoefficient(11) + 5100*poly1_.getCoefficient(10)) + 2286*poly1_.getCoefficient(9)) + 1008*poly1_.getCoefficient(8)) + 434*poly1_.getCoefficient(7)) + 180*poly1_.getCoefficient(6)) + 70*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*(h*(h*(h*(h*(h*(h*(28050*h*poly1_.getCoefficient(11) + 11430*poly1_.getCoefficient(10)) + 4536*poly1_.getCoefficient(9)) + 1736*poly1_.getCoefficient(8)) + 630*poly1_.getCoefficient(7)) + 210*poly1_.getCoefficient(6)) + 60*poly1_.getCoefficient(5)) + 12*poly1_.getCoefficient(4)) + x*(h*h*(h*(h*(h*(h*(h*(41910*h*poly1_.getCoefficient(11) + 15120*poly1_.getCoefficient(10)) + 5208*poly1_.getCoefficient(9)) + 1680*poly1_.getCoefficient(8)) + 490*poly1_.getCoefficient(7)) + 120*poly1_.getCoefficient(6)) + 20*poly1_.getCoefficient(5)) + x*(h*h*(h*(h*(h*(h*(41580*h*poly1_.getCoefficient(11) + 13020*poly1_.getCoefficient(10)) + 3780*poly1_.getCoefficient(9)) + 980*poly1_.getCoefficient(8)) + 210*poly1_.getCoefficient(7)) + 30*poly1_.getCoefficient(6)) + x*(h*h*(h*(h*(h*(28644*h*poly1_.getCoefficient(11) + 7560*poly1_.getCoefficient(10)) + 1764*poly1_.getCoefficient(9)) + 336*poly1_.getCoefficient(8)) + 42*poly1_.getCoefficient(7)) + x*(h*h*(h*(h*(13860*h*poly1_.getCoefficient(11) + 2940*poly1_.getCoefficient(10)) + 504*poly1_.getCoefficient(9)) + 56*poly1_.getCoefficient(8)) + x*(h*h*(h*(4620*h*poly1_.getCoefficient(11) + 720*poly1_.getCoefficient(10)) + 72*poly1_.getCoefficient(9)) + x*(110*h*h*poly1_.getCoefficient(11)*x + h*h*(990*h*poly1_.getCoefficient(11) + 90*poly1_.getCoefficient(10))))))))));
      deltas[3] = h*h*h*(h*(h*(h*(h*(h*(h*(h*(171006*h*poly1_.getCoefficient(11) + 55980*poly1_.getCoefficient(10)) + 18150*poly1_.getCoefficient(9)) + 5796*poly1_.getCoefficient(8)) + 1806*poly1_.getCoefficient(7)) + 540*poly1_.getCoefficient(6)) + 150*poly1_.getCoefficient(5)) + 36*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*h*(h*(h*(h*(h*(h*(h*(615780*h*poly1_.getCoefficient(11) + 181500*poly1_.getCoefficient(10)) + 52164*poly1_.getCoefficient(9)) + 14448*poly1_.getCoefficient(8)) + 3780*poly1_.getCoefficient(7)) + 900*poly1_.getCoefficient(6)) + 180*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(h*h*h*(h*(h*(h*(h*(h*(998250*h*poly1_.getCoefficient(11) + 260820*poly1_.getCoefficient(10)) + 65016*poly1_.getCoefficient(9)) + 15120*poly1_.getCoefficient(8)) + 3150*poly1_.getCoefficient(7)) + 540*poly1_.getCoefficient(6)) + 60*poly1_.getCoefficient(5)) + x*(h*h*h*(h*(h*(h*(h*(956340*h*poly1_.getCoefficient(11) + 216720*poly1_.getCoefficient(10)) + 45360*poly1_.getCoefficient(9)) + 8400*poly1_.getCoefficient(8)) + 1260*poly1_.getCoefficient(7)) + 120*poly1_.getCoefficient(6)) + x*(h*h*h*(h*(h*(h*(595980*h*poly1_.getCoefficient(11) + 113400*poly1_.getCoefficient(10)) + 18900*poly1_.getCoefficient(9)) + 2520*poly1_.getCoefficient(8)) + 210*poly1_.getCoefficient(7)) + x*(h*h*h*(h*(h*(249480*h*poly1_.getCoefficient(11) + 37800*poly1_.getCoefficient(10)) + 4536*poly1_.getCoefficient(9)) + 336*poly1_.getCoefficient(8)) + x*(h*h*h*(h*(69300*h*poly1_.getCoefficient(11) + 7560*poly1_.getCoefficient(10)) + 504*poly1_.getCoefficient(9)) + x*(990*h*h*h*poly1_.getCoefficient(11)*x + h*h*h*(11880*h*poly1_.getCoefficient(11) + 720*poly1_.getCoefficient(10)))))))));
      deltas[4] = h*h*h*h*(h*(h*(h*(h*(h*(h*(3498000*h*poly1_.getCoefficient(11) + 818520*poly1_.getCoefficient(10)) + 186480*poly1_.getCoefficient(9)) + 40824*poly1_.getCoefficient(8)) + 8400*poly1_.getCoefficient(7)) + 1560*poly1_.getCoefficient(6)) + 240*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(h*h*h*h*(h*(h*(h*(h*(h*(9003720*h*poly1_.getCoefficient(11) + 1864800*poly1_.getCoefficient(10)) + 367416*poly1_.getCoefficient(9)) + 67200*poly1_.getCoefficient(8)) + 10920*poly1_.getCoefficient(7)) + 1440*poly1_.getCoefficient(6)) + 120*poly1_.getCoefficient(5)) + x*(h*h*h*h*(h*(h*(h*(h*(10256400*h*poly1_.getCoefficient(11) + 1837080*poly1_.getCoefficient(10)) + 302400*poly1_.getCoefficient(9)) + 43680*poly1_.getCoefficient(8)) + 5040*poly1_.getCoefficient(7)) + 360*poly1_.getCoefficient(6)) + x*(h*h*h*h*(h*(h*(h*(6735960*h*poly1_.getCoefficient(11) + 1008000*poly1_.getCoefficient(10)) + 131040*poly1_.getCoefficient(9)) + 13440*poly1_.getCoefficient(8)) + 840*poly1_.getCoefficient(7)) + x*(h*h*h*h*(h*(h*(2772000*h*poly1_.getCoefficient(11) + 327600*poly1_.getCoefficient(10)) + 30240*poly1_.getCoefficient(9)) + 1680*poly1_.getCoefficient(8)) + x*(h*h*h*h*(h*(720720*h*poly1_.getCoefficient(11) + 60480*poly1_.getCoefficient(10)) + 3024*poly1_.getCoefficient(9)) + x*(7920*h*h*h*h*poly1_.getCoefficient(11)*x + h*h*h*h*(110880*h*poly1_.getCoefficient(11) + 5040*poly1_.getCoefficient(10))))))));
      deltas[5] = h*h*h*h*h*(h*(h*(h*(h*(h*(29607600*h*poly1_.getCoefficient(11) + 5103000*poly1_.getCoefficient(10)) + 834120*poly1_.getCoefficient(9)) + 126000*poly1_.getCoefficient(8)) + 16800*poly1_.getCoefficient(7)) + 1800*poly1_.getCoefficient(6)) + 120*poly1_.getCoefficient(5)) + x*(h*h*h*h*h*(h*(h*(h*(h*(56133000*h*poly1_.getCoefficient(11) + 8341200*poly1_.getCoefficient(10)) + 1134000*poly1_.getCoefficient(9)) + 134400*poly1_.getCoefficient(8)) + 12600*poly1_.getCoefficient(7)) + 720*poly1_.getCoefficient(6)) + x*(h*h*h*h*h*(h*(h*(h*(45876600*h*poly1_.getCoefficient(11) + 5670000*poly1_.getCoefficient(10)) + 604800*poly1_.getCoefficient(9)) + 50400*poly1_.getCoefficient(8)) + 2520*poly1_.getCoefficient(7)) + x*(h*h*h*h*h*(h*(h*(20790000*h*poly1_.getCoefficient(11) + 2016000*poly1_.getCoefficient(10)) + 151200*poly1_.getCoefficient(9)) + 6720*poly1_.getCoefficient(8)) + x*(h*h*h*h*h*(h*(5544000*h*poly1_.getCoefficient(11) + 378000*poly1_.getCoefficient(10)) + 15120*poly1_.getCoefficient(9)) + x*(55440*h*h*h*h*h*poly1_.getCoefficient(11)*x + h*h*h*h*h*(831600*h*poly1_.getCoefficient(11) + 30240*poly1_.getCoefficient(10)))))));
      deltas[6] = h*h*h*h*h*h*(h*(h*(h*(h*(129230640*h*poly1_.getCoefficient(11) + 16435440*poly1_.getCoefficient(10)) + 1905120*poly1_.getCoefficient(9)) + 191520*poly1_.getCoefficient(8)) + 15120*poly1_.getCoefficient(7)) + 720*poly1_.getCoefficient(6)) + x*(h*h*h*h*h*h*(h*(h*(h*(180789840*h*poly1_.getCoefficient(11) + 19051200*poly1_.getCoefficient(10)) + 1723680*poly1_.getCoefficient(9)) + 120960*poly1_.getCoefficient(8)) + 5040*poly1_.getCoefficient(7)) + x*(h*h*h*h*h*h*(h*(h*(104781600*h*poly1_.getCoefficient(11) + 8618400*poly1_.getCoefficient(10)) + 544320*poly1_.getCoefficient(9)) + 20160*poly1_.getCoefficient(8)) + x*(h*h*h*h*h*h*(h*(31600800*h*poly1_.getCoefficient(11) + 1814400*poly1_.getCoefficient(10)) + 60480*poly1_.getCoefficient(9)) + x*(332640*h*h*h*h*h*h*poly1_.getCoefficient(11)*x + h*h*h*h*h*h*(4989600*h*poly1_.getCoefficient(11) + 151200*poly1_.getCoefficient(10))))));
      deltas[7] = h*h*h*h*h*h*h*(h*(h*(h*(322494480*h*poly1_.getCoefficient(11) + 29635200*poly1_.getCoefficient(10)) + 2328480*poly1_.getCoefficient(9)) + 141120*poly1_.getCoefficient(8)) + 5040*poly1_.getCoefficient(7)) + x*(h*h*h*h*h*h*h*(h*(h*(325987200*h*poly1_.getCoefficient(11) + 23284800*poly1_.getCoefficient(10)) + 1270080*poly1_.getCoefficient(9)) + 40320*poly1_.getCoefficient(8)) + x*(h*h*h*h*h*h*h*(h*(128066400*h*poly1_.getCoefficient(11) + 6350400*poly1_.getCoefficient(10)) + 181440*poly1_.getCoefficient(9)) + x*(1663200*h*h*h*h*h*h*h*poly1_.getCoefficient(11)*x + h*h*h*h*h*h*h*(23284800*h*poly1_.getCoefficient(11) + 604800*poly1_.getCoefficient(10)))));
      deltas[8] = h*h*h*h*h*h*h*h*(h*(h*(479001600*h*poly1_.getCoefficient(11) + 30240000*poly1_.getCoefficient(10)) + 1451520*poly1_.getCoefficient(9)) + 40320*poly1_.getCoefficient(8)) + x*(h*h*h*h*h*h*h*h*(h*(332640000*h*poly1_.getCoefficient(11) + 14515200*poly1_.getCoefficient(10)) + 362880*poly1_.getCoefficient(9)) + x*(6652800*h*h*h*h*h*h*h*h*poly1_.getCoefficient(11)*x + h*h*h*h*h*h*h*h*(79833600*h*poly1_.getCoefficient(11) + 1814400*poly1_.getCoefficient(10))));
      deltas[9] = h*h*h*h*h*h*h*h*h*(h*(419126400*h*poly1_.getCoefficient(11) + 16329600*poly1_.getCoefficient(10)) + 362880*poly1_.getCoefficient(9)) + x*(19958400*h*h*h*h*h*h*h*h*h*poly1_.getCoefficient(11)*x + h*h*h*h*h*h*h*h*h*(179625600*h*poly1_.getCoefficient(11) + 3628800*poly1_.getCoefficient(10)));
      deltas[10] = 39916800*h*h*h*h*h*h*h*h*h*h*poly1_.getCoefficient(11)*x + h*h*h*h*h*h*h*h*h*h*(199584000*h*poly1_.getCoefficient(11) + 3628800*poly1_.getCoefficient(10));
      deltas[11] = 39916800*h*h*h*h*h*h*h*h*h*h*h*poly1_.getCoefficient(11);
    }
    if (Degree == 12) {
      deltas[0] = poly1_.getCoefficient(0) + poly1_.getCoefficient(1)*x + poly1_.getCoefficient(10)*x*x*x*x*x*x*x*x*x*x + poly1_.getCoefficient(11)*x*x*x*x*x*x*x*x*x*x*x + poly1_.getCoefficient(12)*x*x*x*x*x*x*x*x*x*x*x*x + poly1_.getCoefficient(2)*x*x + poly1_.getCoefficient(3)*x*x*x + poly1_.getCoefficient(4)*x*x*x*x + poly1_.getCoefficient(5)*x*x*x*x*x + poly1_.getCoefficient(6)*x*x*x*x*x*x + poly1_.getCoefficient(7)*x*x*x*x*x*x*x + poly1_.getCoefficient(8)*x*x*x*x*x*x*x*x + poly1_.getCoefficient(9)*x*x*x*x*x*x*x*x*x;
      deltas[1] = h*(h*(h*(h*(h*(h*(h*(h*(h*(h*(h*(h*poly1_.getCoefficient(12) + poly1_.getCoefficient(11)) + poly1_.getCoefficient(10)) + poly1_.getCoefficient(9)) + poly1_.getCoefficient(8)) + poly1_.getCoefficient(7)) + poly1_.getCoefficient(6)) + poly1_.getCoefficient(5)) + poly1_.getCoefficient(4)) + poly1_.getCoefficient(3)) + poly1_.getCoefficient(2)) + poly1_.getCoefficient(1)) + x*(h*(h*(h*(h*(h*(h*(h*(h*(h*(h*(12*h*poly1_.getCoefficient(12) + 11*poly1_.getCoefficient(11)) + 10*poly1_.getCoefficient(10)) + 9*poly1_.getCoefficient(9)) + 8*poly1_.getCoefficient(8)) + 7*poly1_.getCoefficient(7)) + 6*poly1_.getCoefficient(6)) + 5*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*(h*(h*(h*(h*(h*(h*(h*(h*(66*h*poly1_.getCoefficient(12) + 55*poly1_.getCoefficient(11)) + 45*poly1_.getCoefficient(10)) + 36*poly1_.getCoefficient(9)) + 28*poly1_.getCoefficient(8)) + 21*poly1_.getCoefficient(7)) + 15*poly1_.getCoefficient(6)) + 10*poly1_.getCoefficient(5)) + 6*poly1_.getCoefficient(4)) + 3*poly1_.getCoefficient(3)) + x*(h*(h*(h*(h*(h*(h*(h*(h*(220*h*poly1_.getCoefficient(12) + 165*poly1_.getCoefficient(11)) + 120*poly1_.getCoefficient(10)) + 84*poly1_.getCoefficient(9)) + 56*poly1_.getCoefficient(8)) + 35*poly1_.getCoefficient(7)) + 20*poly1_.getCoefficient(6)) + 10*poly1_.getCoefficient(5)) + 4*poly1_.getCoefficient(4)) + x*(h*(h*(h*(h*(h*(h*(h*(495*h*poly1_.getCoefficient(12) + 330*poly1_.getCoefficient(11)) + 210*poly1_.getCoefficient(10)) + 126*poly1_.getCoefficient(9)) + 70*poly1_.getCoefficient(8)) + 35*poly1_.getCoefficient(7)) + 15*poly1_.getCoefficient(6)) + 5*poly1_.getCoefficient(5)) + x*(h*(h*(h*(h*(h*(h*(792*h*poly1_.getCoefficient(12) + 462*poly1_.getCoefficient(11)) + 252*poly1_.getCoefficient(10)) + 126*poly1_.getCoefficient(9)) + 56*poly1_.getCoefficient(8)) + 21*poly1_.getCoefficient(7)) + 6*poly1_.getCoefficient(6)) + x*(h*(h*(h*(h*(h*(924*h*poly1_.getCoefficient(12) + 462*poly1_.getCoefficient(11)) + 210*poly1_.getCoefficient(10)) + 84*poly1_.getCoefficient(9)) + 28*poly1_.getCoefficient(8)) + 7*poly1_.getCoefficient(7)) + x*(h*(h*(h*(h*(792*h*poly1_.getCoefficient(12) + 330*poly1_.getCoefficient(11)) + 120*poly1_.getCoefficient(10)) + 36*poly1_.getCoefficient(9)) + 8*poly1_.getCoefficient(8)) + x*(h*(h*(h*(495*h*poly1_.getCoefficient(12) + 165*poly1_.getCoefficient(11)) + 45*poly1_.getCoefficient(10)) + 9*poly1_.getCoefficient(9)) + x*(h*(h*(220*h*poly1_.getCoefficient(12) + 55*poly1_.getCoefficient(11)) + 10*poly1_.getCoefficient(10)) + x*(12*h*poly1_.getCoefficient(12)*x + h*(66*h*poly1_.getCoefficient(12) + 11*poly1_.getCoefficient(11))))))))))));
      deltas[2] = h*h*(h*(h*(h*(h*(h*(h*(h*(h*(h*(4094*h*poly1_.getCoefficient(12) + 2046*poly1_.getCoefficient(11)) + 1022*poly1_.getCoefficient(10)) + 510*poly1_.getCoefficient(9)) + 254*poly1_.getCoefficient(8)) + 126*poly1_.getCoefficient(7)) + 62*poly1_.getCoefficient(6)) + 30*poly1_.getCoefficient(5)) + 14*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + 2*poly1_.getCoefficient(2)) + x*(h*h*(h*(h*(h*(h*(h*(h*(h*(h*(24552*h*poly1_.getCoefficient(12) + 11242*poly1_.getCoefficient(11)) + 5100*poly1_.getCoefficient(10)) + 2286*poly1_.getCoefficient(9)) + 1008*poly1_.getCoefficient(8)) + 434*poly1_.getCoefficient(7)) + 180*poly1_.getCoefficient(6)) + 70*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*(h*(h*(h*(h*(h*(h*(h*(67452*h*poly1_.getCoefficient(12) + 28050*poly1_.getCoefficient(11)) + 11430*poly1_.getCoefficient(10)) + 4536*poly1_.getCoefficient(9)) + 1736*poly1_.getCoefficient(8)) + 630*poly1_.getCoefficient(7)) + 210*poly1_.getCoefficient(6)) + 60*poly1_.getCoefficient(5)) + 12*poly1_.getCoefficient(4)) + x*(h*h*(h*(h*(h*(h*(h*(h*(112200*h*poly1_.getCoefficient(12) + 41910*poly1_.getCoefficient(11)) + 15120*poly1_.getCoefficient(10)) + 5208*poly1_.getCoefficient(9)) + 1680*poly1_.getCoefficient(8)) + 490*poly1_.getCoefficient(7)) + 120*poly1_.getCoefficient(6)) + 20*poly1_.getCoefficient(5)) + x*(h*h*(h*(h*(h*(h*(h*(125730*h*poly1_.getCoefficient(12) + 41580*poly1_.getCoefficient(11)) + 13020*poly1_.getCoefficient(10)) + 3780*poly1_.getCoefficient(9)) + 980*poly1_.getCoefficient(8)) + 210*poly1_.getCoefficient(7)) + 30*poly1_.getCoefficient(6)) + x*(h*h*(h*(h*(h*(h*(99792*h*poly1_.getCoefficient(12) + 28644*poly1_.getCoefficient(11)) + 7560*poly1_.getCoefficient(10)) + 1764*poly1_.getCoefficient(9)) + 336*poly1_.getCoefficient(8)) + 42*poly1_.getCoefficient(7)) + x*(h*h*(h*(h*(h*(57288*h*poly1_.getCoefficient(12) + 13860*poly1_.getCoefficient(11)) + 2940*poly1_.getCoefficient(10)) + 504*poly1_.getCoefficient(9)) + 56*poly1_.getCoefficient(8)) + x*(h*h*(h*(h*(23760*h*poly1_.getCoefficient(12) + 4620*poly1_.getCoefficient(11)) + 720*poly1_.getCoefficient(10)) + 72*poly1_.getCoefficient(9)) + x*(h*h*(h*(6930*h*poly1_.getCoefficient(12) + 990*poly1_.getCoefficient(11)) + 90*poly1_.getCoefficient(10)) + x*(132*h*h*poly1_.getCoefficient(12)*x + h*h*(1320*h*poly1_.getCoefficient(12) + 110*poly1_.getCoefficient(11)))))))))));
      deltas[3] = h*h*h*(h*(h*(h*(h*(h*(h*(h*(h*(519156*h*poly1_.getCoefficient(12) + 171006*poly1_.getCoefficient(11)) + 55980*poly1_.getCoefficient(10)) + 18150*poly1_.getCoefficient(9)) + 5796*poly1_.getCoefficient(8)) + 1806*poly1_.getCoefficient(7)) + 540*poly1_.getCoefficient(6)) + 150*poly1_.getCoefficient(5)) + 36*poly1_.getCoefficient(4)) + 6*poly1_.getCoefficient(3)) + x*(h*h*h*(h*(h*(h*(h*(h*(h*(h*(2052072*h*poly1_.getCoefficient(12) + 615780*poly1_.getCoefficient(11)) + 181500*poly1_.getCoefficient(10)) + 52164*poly1_.getCoefficient(9)) + 14448*poly1_.getCoefficient(8)) + 3780*poly1_.getCoefficient(7)) + 900*poly1_.getCoefficient(6)) + 180*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(h*h*h*(h*(h*(h*(h*(h*(h*(3694680*h*poly1_.getCoefficient(12) + 998250*poly1_.getCoefficient(11)) + 260820*poly1_.getCoefficient(10)) + 65016*poly1_.getCoefficient(9)) + 15120*poly1_.getCoefficient(8)) + 3150*poly1_.getCoefficient(7)) + 540*poly1_.getCoefficient(6)) + 60*poly1_.getCoefficient(5)) + x*(h*h*h*(h*(h*(h*(h*(h*(3993000*h*poly1_.getCoefficient(12) + 956340*poly1_.getCoefficient(11)) + 216720*poly1_.getCoefficient(10)) + 45360*poly1_.getCoefficient(9)) + 8400*poly1_.getCoefficient(8)) + 1260*poly1_.getCoefficient(7)) + 120*poly1_.getCoefficient(6)) + x*(h*h*h*(h*(h*(h*(h*(2869020*h*poly1_.getCoefficient(12) + 595980*poly1_.getCoefficient(11)) + 113400*poly1_.getCoefficient(10)) + 18900*poly1_.getCoefficient(9)) + 2520*poly1_.getCoefficient(8)) + 210*poly1_.getCoefficient(7)) + x*(h*h*h*(h*(h*(h*(1430352*h*poly1_.getCoefficient(12) + 249480*poly1_.getCoefficient(11)) + 37800*poly1_.getCoefficient(10)) + 4536*poly1_.getCoefficient(9)) + 336*poly1_.getCoefficient(8)) + x*(h*h*h*(h*(h*(498960*h*poly1_.getCoefficient(12) + 69300*poly1_.getCoefficient(11)) + 7560*poly1_.getCoefficient(10)) + 504*poly1_.getCoefficient(9)) + x*(h*h*h*(h*(118800*h*poly1_.getCoefficient(12) + 11880*poly1_.getCoefficient(11)) + 720*poly1_.getCoefficient(10)) + x*(1320*h*h*h*poly1_.getCoefficient(12)*x + h*h*h*(17820*h*poly1_.getCoefficient(12) + 990*poly1_.getCoefficient(11))))))))));
      deltas[4] = h*h*h*h*(h*(h*(h*(h*(h*(h*(h*(14676024*h*poly1_.getCoefficient(12) + 3498000*poly1_.getCoefficient(11)) + 818520*poly1_.getCoefficient(10)) + 186480*poly1_.getCoefficient(9)) + 40824*poly1_.getCoefficient(8)) + 8400*poly1_.getCoefficient(7)) + 1560*poly1_.getCoefficient(6)) + 240*poly1_.getCoefficient(5)) + 24*poly1_.getCoefficient(4)) + x*(h*h*h*h*(h*(h*(h*(h*(h*(h*(41976000*h*poly1_.getCoefficient(12) + 9003720*poly1_.getCoefficient(11)) + 1864800*poly1_.getCoefficient(10)) + 367416*poly1_.getCoefficient(9)) + 67200*poly1_.getCoefficient(8)) + 10920*poly1_.getCoefficient(7)) + 1440*poly1_.getCoefficient(6)) + 120*poly1_.getCoefficient(5)) + x*(h*h*h*h*(h*(h*(h*(h*(h*(54022320*h*poly1_.getCoefficient(12) + 10256400*poly1_.getCoefficient(11)) + 1837080*poly1_.getCoefficient(10)) + 302400*poly1_.getCoefficient(9)) + 43680*poly1_.getCoefficient(8)) + 5040*poly1_.getCoefficient(7)) + 360*poly1_.getCoefficient(6)) + x*(h*h*h*h*(h*(h*(h*(h*(41025600*h*poly1_.getCoefficient(12) + 6735960*poly1_.getCoefficient(11)) + 1008000*poly1_.getCoefficient(10)) + 131040*poly1_.getCoefficient(9)) + 13440*poly1_.getCoefficient(8)) + 840*poly1_.getCoefficient(7)) + x*(h*h*h*h*(h*(h*(h*(20207880*h*poly1_.getCoefficient(12) + 2772000*poly1_.getCoefficient(11)) + 327600*poly1_.getCoefficient(10)) + 30240*poly1_.getCoefficient(9)) + 1680*poly1_.getCoefficient(8)) + x*(h*h*h*h*(h*(h*(6652800*h*poly1_.getCoefficient(12) + 720720*poly1_.getCoefficient(11)) + 60480*poly1_.getCoefficient(10)) + 3024*poly1_.getCoefficient(9)) + x*(h*h*h*h*(h*(1441440*h*poly1_.getCoefficient(12) + 110880*poly1_.getCoefficient(11)) + 5040*poly1_.getCoefficient(10)) + x*(11880*h*h*h*h*poly1_.getCoefficient(12)*x + h*h*h*h*(190080*h*poly1_.getCoefficient(12) + 7920*poly1_.getCoefficient(11)))))))));
      deltas[5] = h*h*h*h*h*(h*(h*(h*(h*(h*(h*(165528000*h*poly1_.getCoefficient(12) + 29607600*poly1_.getCoefficient(11)) + 5103000*poly1_.getCoefficient(10)) + 834120*poly1_.getCoefficient(9)) + 126000*poly1_.getCoefficient(8)) + 16800*poly1_.getCoefficient(7)) + 1800*poly1_.getCoefficient(6)) + 120*poly1_.getCoefficient(5)) + x*(h*h*h*h*h*(h*(h*(h*(h*(h*(355291200*h*poly1_.getCoefficient(12) + 56133000*poly1_.getCoefficient(11)) + 8341200*poly1_.getCoefficient(10)) + 1134000*poly1_.getCoefficient(9)) + 134400*poly1_.getCoefficient(8)) + 12600*poly1_.getCoefficient(7)) + 720*poly1_.getCoefficient(6)) + x*(h*h*h*h*h*(h*(h*(h*(h*(336798000*h*poly1_.getCoefficient(12) + 45876600*poly1_.getCoefficient(11)) + 5670000*poly1_.getCoefficient(10)) + 604800*poly1_.getCoefficient(9)) + 50400*poly1_.getCoefficient(8)) + 2520*poly1_.getCoefficient(7)) + x*(h*h*h*h*h*(h*(h*(h*(183506400*h*poly1_.getCoefficient(12) + 20790000*poly1_.getCoefficient(11)) + 2016000*poly1_.getCoefficient(10)) + 151200*poly1_.getCoefficient(9)) + 6720*poly1_.getCoefficient(8)) + x*(h*h*h*h*h*(h*(h*(62370000*h*poly1_.getCoefficient(12) + 5544000*poly1_.getCoefficient(11)) + 378000*poly1_.getCoefficient(10)) + 15120*poly1_.getCoefficient(9)) + x*(h*h*h*h*h*(h*(13305600*h*poly1_.getCoefficient(12) + 831600*poly1_.getCoefficient(11)) + 30240*poly1_.getCoefficient(10)) + x*(95040*h*h*h*h*h*poly1_.getCoefficient(12)*x + h*h*h*h*h*(1663200*h*poly1_.getCoefficient(12) + 55440*poly1_.getCoefficient(11))))))));
      deltas[6] = h*h*h*h*h*h*(h*(h*(h*(h*(h*(953029440*h*poly1_.getCoefficient(12) + 129230640*poly1_.getCoefficient(11)) + 16435440*poly1_.getCoefficient(10)) + 1905120*poly1_.getCoefficient(9)) + 191520*poly1_.getCoefficient(8)) + 15120*poly1_.getCoefficient(7)) + 720*poly1_.getCoefficient(6)) + x*(h*h*h*h*h*h*(h*(h*(h*(h*(1550767680*h*poly1_.getCoefficient(12) + 180789840*poly1_.getCoefficient(11)) + 19051200*poly1_.getCoefficient(10)) + 1723680*poly1_.getCoefficient(9)) + 120960*poly1_.getCoefficient(8)) + 5040*poly1_.getCoefficient(7)) + x*(h*h*h*h*h*h*(h*(h*(h*(1084739040*h*poly1_.getCoefficient(12) + 104781600*poly1_.getCoefficient(11)) + 8618400*poly1_.getCoefficient(10)) + 544320*poly1_.getCoefficient(9)) + 20160*poly1_.getCoefficient(8)) + x*(h*h*h*h*h*h*(h*(h*(419126400*h*poly1_.getCoefficient(12) + 31600800*poly1_.getCoefficient(11)) + 1814400*poly1_.getCoefficient(10)) + 60480*poly1_.getCoefficient(9)) + x*(h*h*h*h*h*h*(h*(94802400*h*poly1_.getCoefficient(12) + 4989600*poly1_.getCoefficient(11)) + 151200*poly1_.getCoefficient(10)) + x*(665280*h*h*h*h*h*h*poly1_.getCoefficient(12)*x + h*h*h*h*h*h*(11975040*h*poly1_.getCoefficient(12) + 332640*poly1_.getCoefficient(11)))))));
      deltas[7] = h*h*h*h*h*h*h*(h*(h*(h*(h*(3162075840*h*poly1_.getCoefficient(12) + 322494480*poly1_.getCoefficient(11)) + 29635200*poly1_.getCoefficient(10)) + 2328480*poly1_.getCoefficient(9)) + 141120*poly1_.getCoefficient(8)) + 5040*poly1_.getCoefficient(7)) + x*(h*h*h*h*h*h*h*(h*(h*(h*(3869933760*h*poly1_.getCoefficient(12) + 325987200*poly1_.getCoefficient(11)) + 23284800*poly1_.getCoefficient(10)) + 1270080*poly1_.getCoefficient(9)) + 40320*poly1_.getCoefficient(8)) + x*(h*h*h*h*h*h*h*(h*(h*(1955923200*h*poly1_.getCoefficient(12) + 128066400*poly1_.getCoefficient(11)) + 6350400*poly1_.getCoefficient(10)) + 181440*poly1_.getCoefficient(9)) + x*(h*h*h*h*h*h*h*(h*(512265600*h*poly1_.getCoefficient(12) + 23284800*poly1_.getCoefficient(11)) + 604800*poly1_.getCoefficient(10)) + x*(3991680*h*h*h*h*h*h*h*poly1_.getCoefficient(12)*x + h*h*h*h*h*h*h*(69854400*h*poly1_.getCoefficient(12) + 1663200*poly1_.getCoefficient(11))))));
      deltas[8] = h*h*h*h*h*h*h*h*(h*(h*(h*(6411968640*h*poly1_.getCoefficient(12) + 479001600*poly1_.getCoefficient(11)) + 30240000*poly1_.getCoefficient(10)) + 1451520*poly1_.getCoefficient(9)) + 40320*poly1_.getCoefficient(8)) + x*(h*h*h*h*h*h*h*h*(h*(h*(5748019200*h*poly1_.getCoefficient(12) + 332640000*poly1_.getCoefficient(11)) + 14515200*poly1_.getCoefficient(10)) + 362880*poly1_.getCoefficient(9)) + x*(h*h*h*h*h*h*h*h*(h*(1995840000*h*poly1_.getCoefficient(12) + 79833600*poly1_.getCoefficient(11)) + 1814400*poly1_.getCoefficient(10)) + x*(19958400*h*h*h*h*h*h*h*h*poly1_.getCoefficient(12)*x + h*h*h*h*h*h*h*h*(319334400*h*poly1_.getCoefficient(12) + 6652800*poly1_.getCoefficient(11)))));
      deltas[9] = h*h*h*h*h*h*h*h*h*(h*(h*(8083152000*h*poly1_.getCoefficient(12) + 419126400*poly1_.getCoefficient(11)) + 16329600*poly1_.getCoefficient(10)) + 362880*poly1_.getCoefficient(9)) + x*(h*h*h*h*h*h*h*h*h*(h*(5029516800*h*poly1_.getCoefficient(12) + 179625600*poly1_.getCoefficient(11)) + 3628800*poly1_.getCoefficient(10)) + x*(79833600*h*h*h*h*h*h*h*h*h*poly1_.getCoefficient(12)*x + h*h*h*h*h*h*h*h*h*(1077753600*h*poly1_.getCoefficient(12) + 19958400*poly1_.getCoefficient(11))));
      deltas[10] = h*h*h*h*h*h*h*h*h*h*(h*(6187104000*h*poly1_.getCoefficient(12) + 199584000*poly1_.getCoefficient(11)) + 3628800*poly1_.getCoefficient(10)) + x*(239500800*h*h*h*h*h*h*h*h*h*h*poly1_.getCoefficient(12)*x + h*h*h*h*h*h*h*h*h*h*(2395008000*h*poly1_.getCoefficient(12) + 39916800*poly1_.getCoefficient(11)));
      deltas[11] = 479001600*h*h*h*h*h*h*h*h*h*h*h*poly1_.getCoefficient(12)*x + h*h*h*h*h*h*h*h*h*h*h*(2634508800*h*poly1_.getCoefficient(12) + 39916800*poly1_.getCoefficient(11));
      deltas[12] = 479001600*h*h*h*h*h*h*h*h*h*h*h*h*poly1_.getCoefficient(12);
    }
    return deltas[0];
  }

  template<uint_t Degree>
  inline real_t incrementEval(std::vector<real_t>& deltas) {
    static_assert(Degree <= 12, "Polynomial2DEvaluator not implemented for degree larger than 12");
    if (Degree >= 1) {
      deltas[0] += deltas[1];
    }
    if (Degree >= 2) {
      deltas[1] += deltas[2];
    }
    if (Degree >= 3) {
      deltas[2] += deltas[3];
    }
    if (Degree >= 4) {
      deltas[3] += deltas[4];
    }
    if (Degree >= 5) {
      deltas[4] += deltas[5];
    }
    if (Degree >= 6) {
      deltas[5] += deltas[6];
    }
    if (Degree >= 7) {
      deltas[6] += deltas[7];
    }
    if (Degree >= 8) {
      deltas[7] += deltas[8];
    }
    if (Degree >= 9) {
      deltas[8] += deltas[9];
    }
    if (Degree >= 10) {
      deltas[9] += deltas[10];
    }
    if (Degree >= 11) {
      deltas[10] += deltas[11];
    }
    if (Degree >= 12) {
      deltas[11] += deltas[12];
    }
    return deltas[0];
  }

} // namespace polynomialevaluator

class Polynomial2DEvaluator {
public:

  typedef Polynomial1D<MonomialBasis1D> Polynomial1;
  typedef Polynomial2D<MonomialBasis2D> Polynomial2;

  explicit Polynomial2DEvaluator(uint_t degree)
    : degree_(degree)
    , poly1_(degree)
    , deltas_(degree+1)
  {}

  explicit Polynomial2DEvaluator(const Polynomial2& poly)
    : Polynomial2DEvaluator(poly.getDegree())
  {
    setPolynomial(poly);
  }

  void setPolynomial(const Polynomial2& poly)
  {
    WALBERLA_ASSERT(poly.getDegree() == degree_, "Polynomial degrees don't match!");
    poly2_ = std::make_shared<Polynomial2>(poly);
    // std::cout << "ptr to poly2_ = " << poly2_ << std::endl;
  }

  [[nodiscard]] real_t eval(const Point2D &x) const {
    return poly2_->eval(x);
  }

  void setY(real_t y) {
    // std::cout << "Evaluator2D(q=" << degree_ << ")::setY of p2 with q = "<< poly2_->getDegree() << "\n";

    for (uint_t degree = 0; degree <= degree_; ++degree) {
      poly1_.setCoefficient(degree, 0.0);
    }

    uint_t start = 0;
    real_t y_;

    for (uint_t coeff = 0; coeff <= degree_; ++coeff) {

      uint_t idx = start;
      y_ = real_t(1.0);

      for(uint_t degree = 0; degree <= degree_-coeff; ++degree) {

        poly1_.addToCoefficient(coeff, poly2_->getCoefficient(idx) * y_);

        idx += coeff + degree + 2;
        y_ *= y;
      }

      start += coeff + 1;
    }
  }

  real_t evalX(real_t x) const {
    // return poly1_.eval(x);
    real_t px = poly1_.getCoefficient(0);
    real_t xk = 1.0;
    for (uint_t k = 1; k <= degree_; ++k)
    {
      xk *= x;
      px += poly1_.getCoefficient(k) * xk;
    }
    return px;
  }

  real_t setStartX(real_t x, real_t h) {
    switch (degree_)
    {
    case 0: return polynomialevaluator::setStartX<0>(x, h, poly1_, deltas_ );
    case 1: return polynomialevaluator::setStartX<1>(x, h, poly1_, deltas_ );
    case 2: return polynomialevaluator::setStartX<2>(x, h, poly1_, deltas_ );
    case 3: return polynomialevaluator::setStartX<3>(x, h, poly1_, deltas_ );
    case 4: return polynomialevaluator::setStartX<4>(x, h, poly1_, deltas_ );
    case 5: return polynomialevaluator::setStartX<5>(x, h, poly1_, deltas_ );
    case 6: return polynomialevaluator::setStartX<6>(x, h, poly1_, deltas_ );
    case 7: return polynomialevaluator::setStartX<7>(x, h, poly1_, deltas_ );
    case 8: return polynomialevaluator::setStartX<8>(x, h, poly1_, deltas_ );
    case 9: return polynomialevaluator::setStartX<9>(x, h, poly1_, deltas_ );
    case 10: return polynomialevaluator::setStartX<10>(x, h, poly1_, deltas_ );
    case 11: return polynomialevaluator::setStartX<11>(x, h, poly1_, deltas_ );
    case 12: return polynomialevaluator::setStartX<12>(x, h, poly1_, deltas_ );
    default:return 0;
    }
  }

  real_t incrementEval() {
    switch (degree_)
    {
    case 0: return polynomialevaluator::incrementEval<0>( deltas_ );
    case 1: return polynomialevaluator::incrementEval<1>( deltas_ );
    case 2: return polynomialevaluator::incrementEval<2>( deltas_ );
    case 3: return polynomialevaluator::incrementEval<3>( deltas_ );
    case 4: return polynomialevaluator::incrementEval<4>( deltas_ );
    case 5: return polynomialevaluator::incrementEval<5>( deltas_ );
    case 6: return polynomialevaluator::incrementEval<6>( deltas_ );
    case 7: return polynomialevaluator::incrementEval<7>( deltas_ );
    case 8: return polynomialevaluator::incrementEval<8>( deltas_ );
    case 9: return polynomialevaluator::incrementEval<9>( deltas_ );
    case 10: return polynomialevaluator::incrementEval<10>( deltas_ );
    case 11: return polynomialevaluator::incrementEval<11>( deltas_ );
    case 12: return polynomialevaluator::incrementEval<12>( deltas_ );
    default:return 0;
    }
  }

private:
  std::shared_ptr<Polynomial2> poly2_{};

  uint_t degree_;
  Polynomial1 poly1_;

  std::vector<real_t> deltas_;

};



class Polynomial3DEvaluator{
 public:

  typedef Polynomial1D<MonomialBasis1D> Polynomial1;
  typedef Polynomial2D<MonomialBasis2D> Polynomial2;
  typedef Polynomial3D<MonomialBasis3D> Polynomial3;

  Polynomial3DEvaluator(uint_t degree)
    : degree_(degree)
    , poly_z_(degree)
    , poly_yz_(degree)
    , deltas_(degree+1)
  {
  }

  Polynomial3DEvaluator(const Polynomial3& poly)
    : Polynomial3DEvaluator(poly.getDegree())
  {
    setPolynomial(poly);
  }

  void setPolynomial(const Polynomial3& poly)
  {
    WALBERLA_ASSERT(poly.getDegree() == degree_, "Polynomial degrees don't match!");
    poly_ = std::make_shared<Polynomial3>(poly);
  }

  // restrict polynomial p to z
  void setZ(real_t z) {
    poly_z_.setZero();

    // idx of coefficient c_ijk of 3d polynoial
    uint_t idx_ijk = 0;

    // p(x,y,z) = sum_{ |i+j+k| = 0,...,degree_ } c_ijk * x^i*y^j*z^k
    for (uint_t ijk = 0; ijk <= degree_; ++ijk)
    {
      for (int i = int(ijk); i >= 0; --i)
      {
        // idx of coefficient c_ij of 2d polynoial
        uint_t idx_ij = ijk * (ijk+1) / 2 + ijk - i;
        // z^k
        real_t z_k = real_t(1.0);

        /* compute coefficients of 2d polynomial:
          p|_{z}(x,y) = sum{ |i+j| = 0,...,degree_ } ( sum_{ k = 0,...,degree_ - |i+j| } c_ijk * z^k ) * x^i*y^j
        */
        for (uint_t k = 0; k <= ijk - i; ++k) // note: j = ijk - i - k
        {
          // c_ij += c_ijk * z^k
          poly_z_.addToCoefficient(idx_ij, poly_->getCoefficient(idx_ijk) * z_k);

          idx_ij -= (ijk - k + 1);
          ++idx_ijk;
          z_k *= z;
        }
      }
    }
  }

  // restrict polynomial p|z to y
  void setY(real_t y) {
    poly_yz_.setZero();

    uint_t start = 0;
    real_t y_k;

    for (uint_t coeff = 0; coeff <= degree_; ++coeff) {

      uint_t idx = start;
      y_k = real_t(1.0);

      for(uint_t degree = 0; degree <= degree_-coeff; ++degree) {

        poly_yz_.addToCoefficient(coeff, poly_z_.getCoefficient(idx) * y_k);

        idx += coeff + degree + 2;
        y_k *= y;
      }

      start += coeff + 1;
    }
  }

  // evaluate p|yz(x)
  real_t evalX(real_t x) const {
    // return poly_yz_.eval(x);
    real_t px = poly_yz_.getCoefficient(0);
    real_t xk = 1.0;
    for (uint_t k = 1; k <= degree_; ++k)
    {
      xk *= x;
      px += poly_yz_.getCoefficient(k) * xk;
    }
    return px;
  }

  real_t setStartX(real_t x, real_t h) {
    switch (degree_)
    {
    case 0: return polynomialevaluator::setStartX<0>(x, h, poly_yz_, deltas_ );
    case 1: return polynomialevaluator::setStartX<1>(x, h, poly_yz_, deltas_ );
    case 2: return polynomialevaluator::setStartX<2>(x, h, poly_yz_, deltas_ );
    case 3: return polynomialevaluator::setStartX<3>(x, h, poly_yz_, deltas_ );
    case 4: return polynomialevaluator::setStartX<4>(x, h, poly_yz_, deltas_ );
    case 5: return polynomialevaluator::setStartX<5>(x, h, poly_yz_, deltas_ );
    case 6: return polynomialevaluator::setStartX<6>(x, h, poly_yz_, deltas_ );
    case 7: return polynomialevaluator::setStartX<7>(x, h, poly_yz_, deltas_ );
    case 8: return polynomialevaluator::setStartX<8>(x, h, poly_yz_, deltas_ );
    case 9: return polynomialevaluator::setStartX<9>(x, h, poly_yz_, deltas_ );
    case 10: return polynomialevaluator::setStartX<10>(x, h, poly_yz_, deltas_ );
    case 11: return polynomialevaluator::setStartX<11>(x, h, poly_yz_, deltas_ );
    case 12: return polynomialevaluator::setStartX<12>(x, h, poly_yz_, deltas_ );
    default:return 0;
    }
  }

  real_t incrementEval() {
    switch (degree_)
    {
    case 0: return polynomialevaluator::incrementEval<0>( deltas_ );
    case 1: return polynomialevaluator::incrementEval<1>( deltas_ );
    case 2: return polynomialevaluator::incrementEval<2>( deltas_ );
    case 3: return polynomialevaluator::incrementEval<3>( deltas_ );
    case 4: return polynomialevaluator::incrementEval<4>( deltas_ );
    case 5: return polynomialevaluator::incrementEval<5>( deltas_ );
    case 6: return polynomialevaluator::incrementEval<6>( deltas_ );
    case 7: return polynomialevaluator::incrementEval<7>( deltas_ );
    case 8: return polynomialevaluator::incrementEval<8>( deltas_ );
    case 9: return polynomialevaluator::incrementEval<9>( deltas_ );
    case 10: return polynomialevaluator::incrementEval<10>( deltas_ );
    case 11: return polynomialevaluator::incrementEval<11>( deltas_ );
    case 12: return polynomialevaluator::incrementEval<12>( deltas_ );
    default:return 0;
    }
  }

 private:

   std::shared_ptr<Polynomial3>  poly_{};    // 3d polynomial p

  uint_t              degree_;  // polynomial degree
  Polynomial2         poly_z_;  // p restricted to given z, i.e., p|z
  Polynomial1         poly_yz_;  // p|z restricted to given y, i.e., p|yz

  std::vector<real_t> deltas_;

};

}
