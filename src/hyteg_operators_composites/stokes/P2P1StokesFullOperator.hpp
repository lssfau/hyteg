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

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg_operators_composites/divergence/P2ToP1DivergenceOperator.hpp"
#include "hyteg_operators_composites/gradient/P1ToP2GradientOperator.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesOperatorTemplate.hpp"
#include "hyteg_operators_composites/viscousblock/P2ViscousBlockFullOperator.hpp"
#include "hyteg_operators/operators/full_stokes/P2VectorElementwiseFullStokesP1ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/full_stokes/P2VectorElementwiseEpsilonP0ViscosityIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/full_stokes/P2VectorElementwiseFullStokesP0ViscosityAnnulusMap.hpp"

namespace hyteg {
namespace operatorgeneration {

/// Implements the "full" Stokes operator, i.e., the discrete operator that arises from the discretization of
///
///     - ∇ · ( 2μ ε(u) ) - ( 2/dim ) ∇ ( ∇ · u ) + ∇ p = f,
///       ∇ · u                                         = g.
///
/// where
///
///     ε(w) := (1/2) (∇w + (∇w)ᵀ)
///
/// and μ = μ(x) is a user-supplied finite element function.
///
/// It combines the viscous "full" operator for A with the discrete divergence B and gradient B to correspond to the block
/// matrix
///
///         /        \
///     K = |  A  Bᵀ |
///         |  B  0  |
///         \        /
///
using P2P1StokesFullOperator = detail::P2P1StokesVarViscOperatorTemplate< operatorgeneration::P2ViscousBlockFullOperator,
                                                                          operatorgeneration::P1ToP2GradientOperator,
                                                                          operatorgeneration::P2ToP1DivergenceOperator >;

/// P2P1StokesFullOperator with AnnulusMap blending. See documentation of P2P1StokesFullOperator.
using P2P1StokesFullAnnulusMapOperator =
    detail::P2P1StokesVarViscOperatorTemplate< operatorgeneration::P2ViscousBlockFullAnnulusMapOperator,
                                               operatorgeneration::P1ToP2GradientAnnulusMapOperator,
                                               operatorgeneration::P2ToP1DivergenceAnnulusMapOperator >;

using P2P1StokesFullP0ViscosityAnnulusMapOperator =
    detail::P2P1StokesP0VarViscOperatorTemplate< operatorgeneration::P2VectorElementwiseFullStokesP0ViscosityAnnulusMap,
                                               operatorgeneration::P1ToP2GradientAnnulusMapOperator,
                                               operatorgeneration::P2ToP1DivergenceAnnulusMapOperator >;

/// P2P1StokesFullOperator with IcosahedralShellMap blending. See documentation of P2P1StokesFullOperator.
using P2P1StokesFullIcosahedralShellMapOperator =
    detail::P2P1StokesVarViscOperatorTemplate< operatorgeneration::P2ViscousBlockFullIcosahedralShellMapOperator,
                                               operatorgeneration::P1ToP2GradientIcosahedralShellMapOperator,
                                               operatorgeneration::P2ToP1DivergenceIcosahedralShellMapOperator >;

using P2P1StokesFullP1ViscosityIcosahedralShellMapOperator =
    detail::P2P1StokesP1VarViscOperatorTemplate< operatorgeneration::P2VectorElementwiseFullStokesP1ViscosityIcosahedralShellMap,
                                               operatorgeneration::P1ToP2GradientIcosahedralShellMapOperator,
                                               operatorgeneration::P2ToP1DivergenceIcosahedralShellMapOperator >;

using P2P1StokesFullP0ViscosityIcosahedralShellMapOperator =
    detail::P2P1StokesP0VarViscOperatorTemplate< operatorgeneration::P2VectorElementwiseEpsilonP0ViscosityIcosahedralShellMap,
                                               operatorgeneration::P1ToP2GradientIcosahedralShellMapOperator,
                                               operatorgeneration::P2ToP1DivergenceIcosahedralShellMapOperator >;
} // namespace operatorgeneration
} // namespace hyteg
