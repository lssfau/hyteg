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

#include "add_2D_macroface_edgedof_1_rhsfunction.hpp"
#include "add_2D_macroface_edgedof_2_rhsfunctions.hpp"
#include "add_2D_macroface_edgedof_3_rhsfunctions.hpp"
#include "apply_2D_macroface_edgedof_to_edgedof_add.hpp"
#include "apply_2D_macroface_edgedof_to_edgedof_replace.hpp"
#include "apply_3D_macrocell_edgedof_to_edgedof_add.hpp"
#include "apply_3D_macrocell_edgedof_to_edgedof_replace.hpp"
#include "apply_3D_macroface_one_sided_edgedof_to_edgedof_add.hpp"
#include "apply_3D_macroface_one_sided_edgedof_to_edgedof_add_impl.hpp"
#include "apply_3D_macroface_one_sided_edgedof_to_edgedof_replace.hpp"
#include "apply_3D_macroface_one_sided_edgedof_to_edgedof_replace_impl.hpp"
#include "assign_2D_macroface_edgedof_1_rhsfunction.hpp"
#include "assign_2D_macroface_edgedof_2_rhsfunctions.hpp"
#include "assign_2D_macroface_edgedof_3_rhsfunctions.hpp"
#include "assign_3D_macrocell_edgedof_1_rhsfunction.hpp"
#include "assign_3D_macrocell_edgedof_2_rhsfunctions.hpp"
#include "assign_3D_macrocell_edgedof_3_rhsfunctions.hpp"
#include "communicate_buffered_pack_edgedof_face_to_cell.hpp"
#include "communicate_buffered_pack_edgedof_face_to_cell_impl.hpp"
#include "communicate_buffered_unpack_edgedof_face_to_cell.hpp"
#include "communicate_buffered_unpack_edgedof_face_to_cell_impl.hpp"
#include "communicate_directly_edgedof_cell_to_face_part_1.hpp"
#include "communicate_directly_edgedof_cell_to_face_part_1_impl.hpp"
#include "communicate_directly_edgedof_cell_to_face_part_2.hpp"
#include "communicate_directly_edgedof_cell_to_face_part_2_impl.hpp"
#include "communicate_directly_edgedof_face_to_cell.hpp"
#include "communicate_directly_edgedof_face_to_cell_impl.hpp"