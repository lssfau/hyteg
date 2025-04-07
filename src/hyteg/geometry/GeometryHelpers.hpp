/*
 * Copyright (c) 2022 Nils Kohl, Marcus Mohr.
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

#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/types/PointND.hpp"

namespace hyteg {

/// Find macro-face associated with a given point in 2D
std::tuple< bool, PrimitiveID > findFaceIDForPointIn2D( std::shared_ptr< PrimitiveStorage > storage,
                                                        const Point3D&                      computationalCoords,
                                                        real_t                              searchToleranceRadius );

/// Find macro-cell associated with a given point in 3D
std::tuple< bool, PrimitiveID > findCellIDForPointIn3D( std::shared_ptr< PrimitiveStorage > storage,
                                                        const Point3D&                      computationalCoords,
                                                        real_t                              searchToleranceRadius );

} // namespace hyteg
