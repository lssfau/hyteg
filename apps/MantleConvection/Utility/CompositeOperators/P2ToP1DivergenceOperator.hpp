/*
* Copyright (c) 2017-2025 Nils Kohl, Andreas Burkhart.
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

#include "../../Operators/divergence/P2ToP1ElementwiseDivergenceAnnulusMap_0_0.hpp"
#include "../../Operators/divergence/P2ToP1ElementwiseDivergenceAnnulusMap_0_1.hpp"
#include "../../Operators/divergence/P2ToP1ElementwiseDivergenceIcosahedralShellMap_0_0.hpp"
#include "../../Operators/divergence/P2ToP1ElementwiseDivergenceIcosahedralShellMap_0_1.hpp"
#include "../../Operators/divergence/P2ToP1ElementwiseDivergenceIcosahedralShellMap_0_2.hpp"
#include "../../Operators/divergence/P2ToP1ElementwiseDivergence_0_0.hpp"
#include "../../Operators/divergence/P2ToP1ElementwiseDivergence_0_1.hpp"
#include "../../Operators/divergence/P2ToP1ElementwiseDivergence_0_2.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergenceAnnulusMap_0_0.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergenceAnnulusMap_0_1.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_0.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_1.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_2.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergence_0_0.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergence_0_1.hpp"
#include "../../Operators/grad_rho_rho_divergence/P2ToP1ElementwiseGradRhoRhoDivergence_0_2.hpp"
#include "mixed_operator/VectorToScalarOperator.hpp"

#define FUNC_PREFIX

namespace hyteg {

namespace mcoperators {

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
                                                         mcoperators::P2ToP1ElementwiseDivergence_0_0,
                                                         mcoperators::P2ToP1ElementwiseDivergence_0_1,
                                                         mcoperators::P2ToP1ElementwiseDivergence_0_2 >;

/// P2ToP1DivergenceOperator with AnnulusMap blending. See documentation of P2ToP1DivergenceOperator.
using P2ToP1DivergenceAnnulusMapOperator = VectorToScalarOperator< P2VectorFunction,
                                                                   P1Function,
                                                                   mcoperators::P2ToP1ElementwiseDivergenceAnnulusMap_0_0,
                                                                   mcoperators::P2ToP1ElementwiseDivergenceAnnulusMap_0_1,
                                                                   P2ToP1DummyOperator >;

/// P2ToP1DivergenceOperator with IcosahedralShellMap blending. See documentation of P2ToP1DivergenceOperator.
using P2ToP1DivergenceIcosahedralShellMapOperator =
    VectorToScalarOperator< P2VectorFunction,
                            P1Function,
                            mcoperators::P2ToP1ElementwiseDivergenceIcosahedralShellMap_0_0,
                            mcoperators::P2ToP1ElementwiseDivergenceIcosahedralShellMap_0_1,
                            mcoperators::P2ToP1ElementwiseDivergenceIcosahedralShellMap_0_2 >;

///     Divergence + Rho stokes operator for the compressible case.
///
///     Can be used as a BT Block if we want the grad_rho_rho term to be
///     implicit.
///
///     Weak formulation
///
///         - ∫ ( ∇ · u ) v - ∫ ( ∇rho/rho · u ) v
///

using P2ToP1DivergenceALAOperator = VectorToScalarOperator< P2VectorFunction,
                                                            P1Function,
                                                            mcoperators::P2ToP1ElementwiseGradRhoRhoDivergence_0_0,
                                                            mcoperators::P2ToP1ElementwiseGradRhoRhoDivergence_0_1,
                                                            mcoperators::P2ToP1ElementwiseGradRhoRhoDivergence_0_2 >;

/// P2ToP1DivergenceOperator with AnnulusMap blending. See documentation of P2ToP1DivergenceOperator.
using P2ToP1DivergenceALAAnnulusMapOperator =
    VectorToScalarOperator< P2VectorFunction,
                            P1Function,
                            mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceAnnulusMap_0_0,
                            mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceAnnulusMap_0_1,
                            P2ToP1DummyOperator >;

/// P2ToP1DivergenceOperator with IcosahedralShellMap blending. See documentation of P2ToP1DivergenceOperator.
using P2ToP1DivergenceALAIcosahedralShellMapOperator =
    VectorToScalarOperator< P2VectorFunction,
                            P1Function,
                            mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_0,
                            mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_1,
                            mcoperators::P2ToP1ElementwiseGradRhoRhoDivergenceIcosahedralShellMap_0_2 >;

} // namespace mcoperators
} // namespace hyteg