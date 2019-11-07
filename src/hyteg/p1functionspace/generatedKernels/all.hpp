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

#include "add_2D_macroface_vertexdof_1_rhsfunction.hpp"
#include "add_2D_macroface_vertexdof_2_rhsfunctions.hpp"
#include "add_2D_macroface_vertexdof_3_rhsfunctions.hpp"
#include "add_3D_macrocell_vertexdof_1_rhsfunction.hpp"
#include "add_3D_macrocell_vertexdof_1_rhsfunction_colored.hpp"
#include "apply_2D_macroface_vertexdof_to_vertexdof_add.hpp"
#include "apply_2D_macroface_vertexdof_to_vertexdof_replace.hpp"
#include "apply_3D_macrocell_vertexdof_to_vertexdof_add.hpp"
#include "apply_3D_macrocell_vertexdof_to_vertexdof_add_colored_fused.hpp"
#include "apply_3D_macrocell_vertexdof_to_vertexdof_replace.hpp"
#include "apply_3D_macrocell_vertexdof_to_vertexdof_replace_colored_fused.hpp"
#include "apply_3D_macroface_one_sided_vertexdof_to_vertexdof_add.hpp"
#include "apply_3D_macroface_one_sided_vertexdof_to_vertexdof_replace.hpp"
#include "assign_2D_macroface_vertexdof_1_rhsfunction.hpp"
#include "assign_2D_macroface_vertexdof_2_rhsfunctions.hpp"
#include "assign_2D_macroface_vertexdof_3_rhsfunctions.hpp"
#include "assign_3D_macrocell_vertexdof_1_rhsfunction.hpp"
#include "assign_3D_macrocell_vertexdof_1_rhsfunction_colored.hpp"
#include "assign_3D_macrocell_vertexdof_2_rhsfunctions.hpp"
#include "assign_3D_macrocell_vertexdof_2_rhsfunctions_colored.hpp"
#include "assign_3D_macrocell_vertexdof_3_rhsfunctions.hpp"
#include "assign_3D_macrocell_vertexdof_3_rhsfunctions_colored.hpp"
#include "communicate_directly_vertexdof_cell_to_face.hpp"
#include "communicate_directly_vertexdof_cell_to_face_colored.hpp"
#include "communicate_directly_vertexdof_cell_to_face_colored_impl.hpp"
#include "communicate_directly_vertexdof_cell_to_face_impl.hpp"
#include "communicate_directly_vertexdof_face_to_cell.hpp"
#include "communicate_directly_vertexdof_face_to_cell_colored.hpp"
#include "communicate_directly_vertexdof_face_to_cell_colored_impl.hpp"
#include "communicate_directly_vertexdof_face_to_cell_impl.hpp"
#include "gaussseidel_3D_macrocell_P1.hpp"
#include "gaussseidel_3D_macrocell_P1_colored.hpp"
#include "gaussseidel_3D_macrocell_P1_colored_impl.hpp"
#include "sor_2D_macroface_vertexdof_to_vertexdof.hpp"
#include "sor_2D_macroface_vertexdof_to_vertexdof_backwards.hpp"
#include "sor_3D_macrocell_P1.hpp"
#include "sor_3D_macrocell_P1_backwards.hpp"
#include "sor_3D_macrocell_P1_colored.hpp"
#include "sor_3D_macrocell_P1_colored_impl.hpp"
#include "sor_3D_macroface_P1.hpp"
#include "sor_3D_macroface_P1_backwards.hpp"
#include "sor_3D_macroface_P1_backwards_impl.hpp"
#include "sor_3D_macroface_P1_impl.hpp"
#include "sor_3D_macroface_P1_one_sided.hpp"
#include "sor_3D_macroface_P1_one_sided_backwards.hpp"
#include "sor_3D_macroface_P1_one_sided_backwards_impl.hpp"
#include "sor_3D_macroface_P1_one_sided_impl.hpp"