/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#cmakedefine HYTEG_BUILD_WITH_ADIOS2
#cmakedefine HYTEG_BUILD_WITH_MPFR
#cmakedefine HYTEG_BUILD_WITH_PETSC
#cmakedefine HYTEG_PETSC_BUILT_WITH_HDF5
#cmakedefine HYTEG_BUILD_WITH_TRILINOS
#cmakedefine HYTEG_USE_GENERATED_KERNELS
#cmakedefine HYTEG_TERRANEO_MODULE

#ifdef HYTEG_USE_GENERATED_KERNELS
namespace hyteg {
namespace globalDefines {
constexpr bool useGeneratedKernels = true;
} // namespace globalDefines
} // namespace hyteg
#else
namespace hyteg {
namespace globalDefines {
constexpr bool useGeneratedKernels = false;
} // namespace globalDefines
} // namespace hyteg
#endif

// clang-format off
#define HYTEG_ARCH_ENDIANESS ${CMAKE_CXX_BYTE_ORDER}
// clang-format on

#define RESTRICT WALBERLA_RESTRICT