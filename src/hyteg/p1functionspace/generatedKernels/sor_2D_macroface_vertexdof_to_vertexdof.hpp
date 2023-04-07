/*
 * Copyright (c) 2019 Nils Kohl, Dominik Thoennes.
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

//////////////////////////////////////////////////////////////////////////////
// This file is generated! To fix issues, please fix them in the generator. //
//////////////////////////////////////////////////////////////////////////////

#pragma once
#include "core/Macros.h"
#define RESTRICT WALBERLA_RESTRICT

namespace hyteg {
namespace vertexdof {
namespace macroface {
namespace generated {

void sor_2D_macroface_vertexdof_to_vertexdof(double * RESTRICT _data_p1FaceDst, double * RESTRICT _data_p1FaceRhs, double const * RESTRICT const _data_p1FaceStencil, int level, double relax);
void sor_2D_macroface_vertexdof_to_vertexdof(float * RESTRICT _data_p1FaceDst, float * RESTRICT _data_p1FaceRhs, float const * RESTRICT const _data_p1FaceStencil, int level, float relax);

} // namespace generated
} // namespace macroface
} // namespace vertexdof
} // namespace hyteg