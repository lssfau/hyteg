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
#pragma once

#include "prolongate_2D_macroface_P1_push_additive.hpp"
#include "prolongate_2D_macroface_P2_push_from_edgedofs.hpp"
#include "prolongate_2D_macroface_P2_push_from_vertexdofs.hpp"
#include "prolongate_3D_macrocell_P1_push_additive.hpp"
#include "prolongate_3D_macrocell_P1_push_additive_colored.hpp"
#include "prolongate_3D_macrocell_P2_push_from_edgedofs.hpp"
#include "prolongate_3D_macrocell_P2_push_from_edgedofs_level_0_to_1.hpp"
#include "prolongate_3D_macrocell_P2_push_from_vertexdofs.hpp"
#include "restrict_2D_macroface_P1_pull_additive.hpp"
#include "restrict_2D_macroface_P2_update_edgedofs.hpp"
#include "restrict_2D_macroface_P2_update_vertexdofs.hpp"
#include "restrict_3D_macrocell_P1_pull_additive.hpp"
#include "restrict_3D_macrocell_P1_pull_additive_colored.hpp"
#include "restrict_3D_macrocell_P2_update_edgedofs.hpp"
#include "restrict_3D_macrocell_P2_update_edgedofs_level_1_to_0.hpp"
#include "restrict_3D_macrocell_P2_update_vertexdofs.hpp"