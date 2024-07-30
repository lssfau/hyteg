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
#include <string>

#include "hyteg/mesh/HyTeGMeshDir.hpp"

#include "hyteg/HytegDefinitions.hpp"

namespace hyteg {

std::string prependHyTeGMeshDir( const std::string& partialFilePath )
{
   std::string fullPath{ HYTEG_MESH_DIR };
   fullPath += "/" + partialFilePath;
   return fullPath;
}

std::string prependHyTeGMeshDir( const char* partialFilePath )
{
   std::string fullPath{ HYTEG_MESH_DIR };
   fullPath += "/" + std::string{ partialFilePath };
   return fullPath;
}

} // namespace hyteg
