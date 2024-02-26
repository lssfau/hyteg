/*
* Copyright (c) 2017-2024 Nils Kohl.
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

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg_operators/operators/divergence/P2ToP1ElementwiseDivergence_0_0.hpp"
#include "hyteg_operators/operators/divergence/P2ToP1ElementwiseDivergence_0_1.hpp"
#include "hyteg_operators/operators/divergence/P2ToP1ElementwiseDivergence_0_2.hpp"

#include "mixed_operator/VectorToScalarOperator.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

/// Implements the divergence operator from the vectorial P2 space to the scalar P1 space.
/// Specifically, it corresponds to the operator (strong form)
///
///     ∇ · u
///
/// and the weak form (trial: vector valued u, test: scalar v)
///
///     - ∫ ( ∇ · u ) v.
///
/// This gives a 1x3 block operator
///
///         /                \
///     B = | B_11 B_12 B_13 |
///         \                /
///
using P2ToP1DivergenceOperator = VectorToScalarOperator< P2VectorFunction,
                                                         P1Function,
                                                         operatorgeneration::P2ToP1ElementwiseDivergence_0_0,
                                                         operatorgeneration::P2ToP1ElementwiseDivergence_0_1,
                                                         operatorgeneration::P2ToP1ElementwiseDivergence_0_2 >;

} // namespace operatorgeneration
} // namespace hyteg