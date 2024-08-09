/*
 * Copyright (c) 2024 Marcus Mohr.
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

#include "hyteg/HytegDefinitions.hpp"

namespace hyteg {

/// Prepend full path to directory for meshes shipped with HyTeG
///
/// This method takes as input a filename or partial path to a file
/// and prepends to it the full filesystem path to the directory containing
/// the meshes shipped with HyTeG w.r.t. to the build directory.
std::string prependHyTeGMeshDir( const std::string& partialFilePath );

/// Prepend full path to directory for meshes shipped with HyTeG
///
/// This method takes as input a filename or partial path to a file
/// and prepends to it the full filesystem path to the directory containing
/// the meshes shipped with HyTeG w.r.t. to the build directory.
std::string prependHyTeGMeshDir( const char* partialFilePath );

} // namespace hyteg
