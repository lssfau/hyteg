/*
 * Copyright (c) 2017-2021 Nils Kohl.
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

#include "hyteg/geometry/GeometryMap.hpp"
#include "hyteg/forms/form_hyteg_base/P1FormHyTeG.hpp"

namespace hyteg {
namespace forms {
    void MeanAverageCoefficients(std::vector<std::reference_wrapper<real_t>> coefficients) {
        real_t nCoeffs = coefficients.size();
        real_t sum = 0;
        for (auto c : coefficients) {
            sum += c;
        }
        for (auto c : coefficients) {
            c.get() = sum / nCoeffs;
        }
    }
    void HarmonicAverageCoefficients(std::vector<std::reference_wrapper<real_t>> coefficients) {
        real_t nCoeffs = coefficients.size();
        real_t sum = 0;
        for (auto c : coefficients) {
            sum += 1.0/c;
        }
        for (auto c : coefficients) {
            c.get() = 1.0/(sum / nCoeffs);
        }
    }
} // namespace forms
} // namespace hyteg
