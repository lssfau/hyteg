/*
 * Copyright (c) 2017-2024 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#cmakedefine HYTEG_BUILD_WITH_MPI
#cmakedefine HYTEG_BUILD_WITH_ADIOS2
#cmakedefine HYTEG_BUILD_WITH_MPFR
#cmakedefine HYTEG_BUILD_WITH_PETSC
#cmakedefine HYTEG_PETSC_BUILT_WITH_HDF5
#cmakedefine HYTEG_BUILD_WITH_TRILINOS
#cmakedefine HYTEG_USE_GENERATED_KERNELS
#cmakedefine HYTEG_TERRANEO_MODULE
#cmakedefine HYTEG_MANTLECONVECTION_APP
#cmakedefine HYTEG_DONT_BUILD_DURING_CI
#cmakedefine HYTEG_BUILD_WITH_PYTHON3
#cmakedefine HYTEG_USE_SIGNED_INT_FOR_ADIOS2

namespace hyteg {
namespace globalDefines {
#ifdef HYTEG_USE_GENERATED_KERNELS
constexpr bool useGeneratedKernels = true;
#else
constexpr bool useGeneratedKernels = false;
#endif

} // namespace globalDefines
} // namespace hyteg

// clang-format off
#define HYTEG_ARCH_ENDIANESS ${CMAKE_CXX_BYTE_ORDER}
// clang-format on

#ifdef HYTEG_BUILD_WITH_ADIOS2

// set integer type for connectivity information to be used when
// writing BP files for ParaView to either uint64_t or int64_t
#ifdef HYTEG_USE_SIGNED_INT_FOR_ADIOS2
#define ADIOS2_PARAVIEW_INT_TYPE int64_t
#else
#define ADIOS2_PARAVIEW_INT_TYPE uint64_t
#endif

#cmakedefine ADIOS2_CHECKPOINT_FORMAT

#endif

/// Absolute path to the mesh-files shipped with HyTeG
/// in the built-out directory
#define HYTEG_MESH_DIR "${hyteg_BINARY_DIR}/data/meshes"

#define RESTRICT WALBERLA_RESTRICT
