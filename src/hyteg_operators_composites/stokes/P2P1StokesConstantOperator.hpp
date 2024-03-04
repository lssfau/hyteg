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
#include "hyteg_operators_composites/viscousblock/P2ViscousBlockLaplaceOperator.hpp"

namespace hyteg {
namespace operatorgeneration {

/// Implements the constant Stokes operator, i.e., the discrete operator that arises from the discretization of
///
///     - Δ u + ∇ p = f,
///       ∇ · u     = g.
///
/// It combines the block-Laplace operator for A with the discrete divergence B and gradient B to correspond to the block matrix
///
///         /        \
///     K = |  A  Bᵀ |
///         |  B  0  |
///         \        /
///
using P2P1StokesConstantOperator = detail::P2P1StokesConstViscOperatorTemplate< operatorgeneration::P2ViscousBlockLaplaceOperator,
                                                                                operatorgeneration::P1ToP2GradientOperator,
                                                                                operatorgeneration::P2ToP1DivergenceOperator >;

/// P2P1StokesConstantOperator with IcosahedralShellMap blending. See documentation of P2P1StokesConstantOperator.
using P2P1StokesConstantIcosahedralShellMapOperator =
    detail::P2P1StokesConstViscOperatorTemplate< operatorgeneration::P2ViscousBlockLaplaceIcosahedralShellMapOperator,
                                                 operatorgeneration::P1ToP2GradientIcosahedralShellMapOperator,
                                                 operatorgeneration::P2ToP1DivergenceIcosahedralShellMapOperator >;

} // namespace operatorgeneration
} // namespace hyteg
