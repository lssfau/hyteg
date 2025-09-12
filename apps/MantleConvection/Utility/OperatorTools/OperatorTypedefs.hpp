/*
 * Copyright (c) 2024-2025 Andreas Burkhart.
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

#include "../LHS/AdvectionDiffusionOperator.hpp"
#include "../LHS/SaddlePointOperator.hpp"
#include "../RHS/AdvectionDiffusionRHS.hpp"
#include "../RHS/SaddlePointOperatorRHS.hpp"
#include "OptimisationTypedefs.hpp"

// #########################################
// #### basic TALA with frozen velocity ####
// #########################################

// LHS
typedef MantleConvection::SaddlePointOperator< MC_ABlock_Vec, MC_BTBlock, MC_BBlock, MC_Projection, MC_NoOp >
    stokesFrozenVelocityTALA_LHS;

typedef MantleConvection::
    SaddlePointOperator< MC_ABlock_Vec_AnnulusMap, MC_BTBlock_AnnulusMap, MC_BBlock_AnnulusMap, MC_Projection, MC_NoOp >
        stokesFrozenVelocityTALA_AnnulusMap_LHS;

typedef MantleConvection::SaddlePointOperator< MC_ABlock_Vec_IcosahedralShellMap,
                                               MC_BTBlock_IcosahedralShellMap,
                                               MC_BBlock_IcosahedralShellMap,
                                               MC_Projection,
                                               MC_NoOp >
    stokesFrozenVelocityTALA_IcosahedralShellMap_LHS;

// RHS

typedef MantleConvection::SaddlePointOperatorRHS< MC_NoOp,
                                                  MC_NoOp,
                                                  MC_TemperatureToVelocityRHS,
                                                  MC_VelocityToPressureRHS,
                                                  MC_NoOp,
                                                  MC_NoOp,
                                                  MC_NoOp,
                                                  MC_Projection >
    stokesFrozenVelocityTALA_RHS;

typedef MantleConvection::SaddlePointOperatorRHS< MC_NoOp,
                                                  MC_NoOp,
                                                  MC_TemperatureToVelocityRHS_AnnulusMap,
                                                  MC_VelocityToPressureRHS_AnnulusMap,
                                                  MC_NoOp,
                                                  MC_NoOp,
                                                  MC_NoOp,
                                                  MC_Projection >
    stokesFrozenVelocityTALA_AnnulusMap_RHS;

typedef MantleConvection::SaddlePointOperatorRHS< MC_NoOp,
                                                  MC_NoOp,
                                                  MC_TemperatureToVelocityRHS_IcosahedralShellMap,
                                                  MC_VelocityToPressureRHS_IcosahedralShellMap,
                                                  MC_NoOp,
                                                  MC_NoOp,
                                                  MC_NoOp,
                                                  MC_Projection >
    stokesFrozenVelocityTALA_IcosahedralShellMap_RHS;

// ###########################################
// #### SUPG Advection Diffusion Operator ####
// ###########################################

// LHS

typedef MantleConvection::AdvectionDiffusionOperator< MC_P2Mass,
                                                      MC_Advection,
                                                      MC_DivKGrad,
                                                      MC_DiffusionAdditional,
                                                      MC_AdiabaticHeating,

                                                      MC_P2MassSUPG,
                                                      MC_AdvectionSUPG,
                                                      MC_DiffusionSUPG,
                                                      MC_AdiabaticHeatingSUPG >
    transportSUPG_P2_LHS;

typedef MantleConvection::AdvectionDiffusionOperator< MC_P2Mass_AnnulusMap,
                                                      MC_Advection_AnnulusMap,
                                                      MC_DivKGrad_AnnulusMap,
                                                      MC_DiffusionAdditional_AnnulusMap,
                                                      MC_AdiabaticHeating_AnnulusMap,

                                                      MC_P2MassSUPG_AnnulusMap,
                                                      MC_AdvectionSUPG_AnnulusMap,
                                                      MC_DiffusionSUPG_AnnulusMap,
                                                      MC_AdiabaticHeatingSUPG_AnnulusMap >
    transportSUPG_P2_AnnulusMap_LHS;

typedef MantleConvection::AdvectionDiffusionOperator< MC_P2Mass_IcosahedralShellMap,
                                                      MC_Advection_IcosahedralShellMap,
                                                      MC_DivKGrad_IcosahedralShellMap,
                                                      MC_DiffusionAdditional_IcosahedralShellMap,
                                                      MC_AdiabaticHeating_IcosahedralShellMap,

                                                      MC_P2MassSUPG_IcosahedralShellMap,
                                                      MC_AdvectionSUPG_IcosahedralShellMap,
                                                      MC_DiffusionSUPG_IcosahedralShellMap,
                                                      MC_AdiabaticHeatingSUPG_IcosahedralShellMap >
    transportSUPG_P2_IcosahedralShellMap_LHS;

// RHS

typedef MantleConvection::AdvectionDiffusionOperatorRHS< MC_P2Mass,
                                                         MC_P2Mass,
                                                         MC_ShearHeating,
                                                         MC_NoOp,

                                                         MC_P2MassSUPG,
                                                         MC_P2MassSUPG,
                                                         MC_ShearHeatingSUPG,
                                                         MC_NoOp >
    transportSUPG_P2_RHS;

typedef MantleConvection::AdvectionDiffusionOperatorRHS< MC_P2Mass_AnnulusMap,
                                                         MC_P2Mass_AnnulusMap,
                                                         MC_ShearHeating_AnnulusMap,
                                                         MC_NoOp,

                                                         MC_P2MassSUPG_AnnulusMap,
                                                         MC_P2MassSUPG_AnnulusMap,
                                                         MC_ShearHeatingSUPG_AnnulusMap,
                                                         MC_NoOp >
    transportSUPG_P2_AnnulusMap_RHS;

typedef MantleConvection::AdvectionDiffusionOperatorRHS< MC_P2Mass_IcosahedralShellMap,
                                                         MC_P2Mass_IcosahedralShellMap,
                                                         MC_ShearHeating_IcosahedralShellMap,
                                                         MC_NoOp,

                                                         MC_P2MassSUPG_IcosahedralShellMap,
                                                         MC_P2MassSUPG_IcosahedralShellMap,
                                                         MC_ShearHeatingSUPG_IcosahedralShellMap,
                                                         MC_NoOp >
    transportSUPG_P2_IcosahedralShellMap_RHS;

// ###########################################
// #### MMOC Advection Diffusion Operator ####
// ###########################################

// LHS

typedef MantleConvection::
    AdvectionDiffusionOperator< MC_P2Mass, MC_NoOp, MC_DivKGrad, MC_DiffusionAdditional, MC_AdiabaticHeating >
        transportMMOC_P2_LHS;

typedef MantleConvection::AdvectionDiffusionOperator< MC_P2Mass_AnnulusMap,
                                                      MC_NoOp,
                                                      MC_DivKGrad_AnnulusMap,
                                                      MC_DiffusionAdditional_AnnulusMap,
                                                      MC_AdiabaticHeating_AnnulusMap >
    transportMMOC_P2_AnnulusMap_LHS;

typedef MantleConvection::AdvectionDiffusionOperator< MC_P2Mass_IcosahedralShellMap,
                                                      MC_NoOp,
                                                      MC_DivKGrad_IcosahedralShellMap,
                                                      MC_DiffusionAdditional_IcosahedralShellMap,
                                                      MC_AdiabaticHeating_IcosahedralShellMap >
    transportMMOC_P2_IcosahedralShellMap_LHS;

// RHS

typedef MantleConvection::AdvectionDiffusionOperatorRHS< MC_P2Mass, MC_P2Mass, MC_ShearHeating, MC_NoOp > transportMMOC_P2_RHS;

typedef MantleConvection::
    AdvectionDiffusionOperatorRHS< MC_P2Mass_AnnulusMap, MC_P2Mass_AnnulusMap, MC_ShearHeating_AnnulusMap, MC_NoOp >
        transportMMOC_P2_AnnulusMap_RHS;

typedef MantleConvection::AdvectionDiffusionOperatorRHS< MC_P2Mass_IcosahedralShellMap,
                                                         MC_P2Mass_IcosahedralShellMap,
                                                         MC_ShearHeating_IcosahedralShellMap,
                                                         MC_NoOp >
    transportMMOC_P2_IcosahedralShellMap_RHS;