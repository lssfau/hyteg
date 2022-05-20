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

#include "PETScWrapper.hpp"

#include "PETScLUSolver.hpp"

#ifdef HYTEG_BUILD_WITH_PETSC

namespace hyteg {

template<typename ValueType, template <class> class F, class O>
class PETScPreconditioner {
public:
  PETScPreconditioner( O& opr, std::shared_ptr< F< idx_t > >& numerator, uint_t localSize, uint_t globalSize )
  : opr_( opr )
  , petscSolver( numerator, localSize, globalSize )
  {}

  // y = M^{-1} * x
  void apply(F<ValueType> &x, F<ValueType> &y, uint_t level, DoFType flag) {
    petscSolver.solve(opr_,y,x,x,level,0,0,flag);
  }

private:
  O opr_;
  PETScLUSolver<ValueType, F, O> petscSolver;
};

}

#endif