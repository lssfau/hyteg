/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

namespace hyteg {

template<class F, class O1, class O2 >
class StokesBlockDiagonalApplyPreconditioner {
public:

    StokesBlockDiagonalApplyPreconditioner(O1 velocityOpr, O2 pressureOpr, const std::shared_ptr<PrimitiveStorage> & storage, size_t minLevel, size_t maxLevel)
      : velocityOpr_(velocityOpr), pressureOpr_(pressureOpr),
        r("r", storage, minLevel, maxLevel)
  {}

  // y = M^{-1} * x
  void apply(F &x, F &y, uint_t level, DoFType flag)
  {

    y.assign({1.0}, {&x}, level, flag);

    velocityOpr_.apply( x.u, y.u, level, flag, Replace );
    velocityOpr_.apply( x.v, y.v, level, flag, Replace );
    pressureOpr_.apply( x.p, y.p, level, flag | DirichletBoundary, Replace );


  }

private:
  O1 velocityOpr_;
  O2 pressureOpr_;

  F r;
};

}