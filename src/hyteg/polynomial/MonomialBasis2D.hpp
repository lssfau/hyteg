/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes.
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

// This file was generated by the monomial_basis_2d.py Python script
// Do not edit it by hand
#include "hyteg/types/PointND.hpp"

#include "core/DataTypes.h"
#include "core/Abort.h"

namespace hyteg {

class MonomialBasis2D {
public:
  static walberla::real_t eval(walberla::uint_t basis, const Point2D &x) {
    switch(basis) {
      case 0:
        return 1;
      case 1:
        return x[0];
      case 2:
        return x[1];
      case 3:
        return x[0]*x[0];
      case 4:
        return x[0]*x[1];
      case 5:
        return x[1]*x[1];
      case 6:
        return x[0]*x[0]*x[0];
      case 7:
        return x[1]*(x[0]*x[0]);
      case 8:
        return x[0]*(x[1]*x[1]);
      case 9:
        return x[1]*x[1]*x[1];
      case 10:
        return x[0]*x[0]*x[0]*x[0];
      case 11:
        return x[1]*(x[0]*x[0]*x[0]);
      case 12:
        return (x[0]*x[0])*(x[1]*x[1]);
      case 13:
        return x[0]*(x[1]*x[1]*x[1]);
      case 14:
        return x[1]*x[1]*x[1]*x[1];
      case 15:
        return x[0]*x[0]*x[0]*x[0]*x[0];
      case 16:
        return x[1]*(x[0]*x[0]*x[0]*x[0]);
      case 17:
        return (x[1]*x[1])*(x[0]*x[0]*x[0]);
      case 18:
        return (x[0]*x[0])*(x[1]*x[1]*x[1]);
      case 19:
        return x[0]*(x[1]*x[1]*x[1]*x[1]);
      case 20:
        return x[1]*x[1]*x[1]*x[1]*x[1];
      case 21:
        return x[0]*x[0]*x[0]*x[0]*x[0]*x[0];
      case 22:
        return x[1]*(x[0]*x[0]*x[0]*x[0]*x[0]);
      case 23:
        return (x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]);
      case 24:
        return (x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]);
      case 25:
        return (x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]);
      case 26:
        return x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]);
      case 27:
        return x[1]*x[1]*x[1]*x[1]*x[1]*x[1];
      case 28:
        return x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0];
      case 29:
        return x[1]*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 30:
        return (x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]);
      case 31:
        return (x[1]*x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]);
      case 32:
        return (x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]);
      case 33:
        return (x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]);
      case 34:
        return x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 35:
        return x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1];
      case 36:
        return x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0];
      case 37:
        return x[1]*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 38:
        return (x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 39:
        return (x[1]*x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]);
      case 40:
        return (x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]);
      case 41:
        return (x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]);
      case 42:
        return (x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 43:
        return x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 44:
        return x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1];
      case 45:
        return x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0];
      case 46:
        return x[1]*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 47:
        return (x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 48:
        return (x[1]*x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 49:
        return (x[1]*x[1]*x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]);
      case 50:
        return (x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]);
      case 51:
        return (x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 52:
        return (x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 53:
        return x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 54:
        return x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1];
      case 55:
        return x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0];
      case 56:
        return x[1]*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 57:
        return (x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 58:
        return (x[1]*x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 59:
        return (x[1]*x[1]*x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 60:
        return (x[0]*x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]);
      case 61:
        return (x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 62:
        return (x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 63:
        return (x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 64:
        return x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 65:
        return x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1];
      case 66:
        return x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0];
      case 67:
        return x[1]*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 68:
        return (x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 69:
        return (x[1]*x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 70:
        return (x[1]*x[1]*x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 71:
        return (x[1]*x[1]*x[1]*x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 72:
        return (x[0]*x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 73:
        return (x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 74:
        return (x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 75:
        return (x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 76:
        return x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 77:
        return x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1];
      case 78:
        return x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0];
      case 79:
        return x[1]*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 80:
        return (x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 81:
        return (x[1]*x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 82:
        return (x[1]*x[1]*x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 83:
        return (x[1]*x[1]*x[1]*x[1]*x[1])*(x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]);
      case 84:
        return (x[0]*x[0]*x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 85:
        return (x[0]*x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 86:
        return (x[0]*x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 87:
        return (x[0]*x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 88:
        return (x[0]*x[0])*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 89:
        return x[0]*(x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]);
      case 90:
        return x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1]*x[1];
      default:
      WALBERLA_ABORT("Polynomial basis " << basis << " was not generated");
    }
  }
};

}