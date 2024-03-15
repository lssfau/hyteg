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

#include "hyteg/mixedoperators/MixedDummyOperators.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg_operators/operators/gradient/P1ToP2ElementwiseGradientAnnulusMap_0_0.hpp"
#include "hyteg_operators/operators/gradient/P1ToP2ElementwiseGradientAnnulusMap_1_0.hpp"
#include "hyteg_operators/operators/gradient/P1ToP2ElementwiseGradientIcosahedralShellMap_0_0.hpp"
#include "hyteg_operators/operators/gradient/P1ToP2ElementwiseGradientIcosahedralShellMap_1_0.hpp"
#include "hyteg_operators/operators/gradient/P1ToP2ElementwiseGradientIcosahedralShellMap_2_0.hpp"
#include "hyteg_operators/operators/gradient/P1ToP2ElementwiseGradient_0_0.hpp"
#include "hyteg_operators/operators/gradient/P1ToP2ElementwiseGradient_1_0.hpp"
#include "hyteg_operators/operators/gradient/P1ToP2ElementwiseGradient_2_0.hpp"

#include "mixed_operator/ScalarToVectorOperator.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace operatorgeneration {

/// Implements the gradient operator from the scalar P1 space to the vectorial P2 space.
/// Specifically, it corresponds to the operator (strong form)
///
///     ∇u
///
/// and the weak form (trial: scalar u, test: vector valued v)
///
///     - ∫ ( ∇ · v ) u.
///
/// This gives a 3x1 block operator
///
///         /      \
///         | B_11 |
///     B = | B_21 |
///         | B_31 |
///         \      /
///
using P1ToP2GradientOperator = ScalarToVectorOperator< P1Function,
                                                       P2VectorFunction,
                                                       operatorgeneration::P1ToP2ElementwiseGradient_0_0,
                                                       operatorgeneration::P1ToP2ElementwiseGradient_1_0,
                                                       operatorgeneration::P1ToP2ElementwiseGradient_2_0 >;

/// P1ToP2GradientOperator with AnnulusMap blending. See documentation of P1ToP2GradientOperator.
using P1ToP2GradientAnnulusMapOperator = ScalarToVectorOperator< P1Function,
                                                                 P2VectorFunction,
                                                                 operatorgeneration::P1ToP2ElementwiseGradientAnnulusMap_0_0,
                                                                 operatorgeneration::P1ToP2ElementwiseGradientAnnulusMap_1_0,
                                                                 P1ToP2DummyOperator >;

/// P1ToP2GradientOperator with IcosahedralShellMap blending. See documentation of P1ToP2GradientOperator.
using P1ToP2GradientIcosahedralShellMapOperator =
    ScalarToVectorOperator< P1Function,
                            P2VectorFunction,
                            operatorgeneration::P1ToP2ElementwiseGradientIcosahedralShellMap_0_0,
                            operatorgeneration::P1ToP2ElementwiseGradientIcosahedralShellMap_1_0,
                            operatorgeneration::P1ToP2ElementwiseGradientIcosahedralShellMap_2_0 >;

} // namespace operatorgeneration
} // namespace hyteg