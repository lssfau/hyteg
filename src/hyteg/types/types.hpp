/*
 * Copyright (c) 2017-2019 Boerge Struempfel, Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include <cstddef>
#include <iostream>

namespace hyteg {

using idx_t = int64_t;

enum UpdateType
{
  Replace = 0,
  Add = 1
};

enum DoFType:size_t
{
  None = 0,
  All = 1+2+4+8,
  Boundary = 2+4+8,
  Inner = 1,
  DirichletBoundary = 2,
  NeumannBoundary = 4,
  FreeslipBoundary = 8
};

inline DoFType operator|(DoFType a, DoFType b){
  return DoFType(static_cast<size_t>(a) | static_cast<size_t>(b));
}

inline DoFType operator&(DoFType a, DoFType b){
  return DoFType(static_cast<size_t>(a) & static_cast<size_t>(b));
}

inline DoFType operator^(DoFType a, DoFType b){
  return DoFType(static_cast<size_t>(a) ^ static_cast<size_t>(b));
}

inline std::ostream& operator<<(std::ostream &os, const DoFType type){
  switch(type) {
    case NeumannBoundary   :
      return os << "NeumannBoundary";
    case DirichletBoundary :
      return os << "DirichletBoundary";
    default:
      return os << "Inner";
  }
}


inline bool testFlag(DoFType a, DoFType b)
{
  return (a & b) != 0;
}

enum class MemoryType {
  Base,
  P1Stencil,
  P1Function,
  P1BubbleFunction,
  P1BubbleStencil
};

enum class CycleType
{
   VCYCLE,
   WCYCLE
};

} // namespace hyteg

