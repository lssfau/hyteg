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

#include "hyteg/operators/NoOperator.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg_operators/operators/adiabatic_heating/P2ElementwiseAdiabaticHeatingAnnulusMap.hpp"
#include "hyteg_operators/operators/adiabatic_heating/P2ElementwiseAdiabaticHeatingIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/advection/P2ElementwiseAdvection.hpp"
#include "hyteg_operators/operators/advection/P2ElementwiseAdvectionAnnulusMap.hpp"
#include "hyteg_operators/operators/advection/P2ElementwiseAdvectionIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/diffusion_inv_rho/P2ElementwiseDiffusionInvRho.hpp"
#include "hyteg_operators/operators/diffusion_inv_rho/P2ElementwiseDiffusionInvRhoAnnulusMap.hpp"
#include "hyteg_operators/operators/diffusion_inv_rho/P2ElementwiseDiffusionInvRhoIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/div_k_grad/P1ElementwiseDivKGrad.hpp"
#include "hyteg_operators/operators/div_k_grad/P1ElementwiseDivKGradAnnulusMap.hpp"
#include "hyteg_operators/operators/div_k_grad/P1ElementwiseDivKGradIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGradP1Coefficient.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGradP1CoefficientAnnulusMap.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGradP1CoefficientIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/supg_adiabatic_heating/P2ElementwiseSupgAdiabaticHeatingAnnulusMap.hpp"
#include "hyteg_operators/operators/supg_adiabatic_heating/P2ElementwiseSupgAdiabaticHeatingIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/supg_advection/P2ElementwiseSupgAdvection.hpp"
#include "hyteg_operators/operators/supg_advection/P2ElementwiseSupgAdvectionAnnulusMap.hpp"
#include "hyteg_operators/operators/supg_advection/P2ElementwiseSupgAdvectionIcosahedralShellMap.hpp"
#include "hyteg_operators/operators/supg_diffusion/P2ElementwiseSupgDiffusion.hpp"
#include "hyteg_operators/operators/supg_diffusion/P2ElementwiseSupgDiffusionAnnulusMap.hpp"
#include "hyteg_operators/operators/supg_diffusion/P2ElementwiseSupgDiffusionIcosahedralShellMap.hpp"

#include "../../Operators/full_stokes_scaled_divdiv_vectorial/P2VectorElementwiseFullStokesScaledDivdivVectorial.hpp"
#include "../../Operators/full_stokes_scaled_divdiv_vectorial/P2VectorElementwiseFullStokesScaledDivdivVectorialAnnulusMap.hpp"
#include "../../Operators/full_stokes_scaled_divdiv_vectorial/P2VectorElementwiseFullStokesScaledDivdivVectorialIcosahedralShellMap.hpp"
#include "../../Operators/full_stokes_scaled_divdiv_vectorial_frank_kamenetskii_simple_visc_base/P2VectorElementwiseFullStokesScaledDivdivVectorialFrankKamenetskiiSimpleViscBase.hpp"
#include "../../Operators/full_stokes_scaled_divdiv_vectorial_frank_kamenetskii_simple_visc_base/P2VectorElementwiseFullStokesScaledDivdivVectorialFrankKamenetskiiSimpleViscBaseAnnulusMap.hpp"
#include "../../Operators/full_stokes_scaled_divdiv_vectorial_frank_kamenetskii_simple_visc_base/P2VectorElementwiseFullStokesScaledDivdivVectorialFrankKamenetskiiSimpleViscBaseIcosahedralShellMap.hpp"
#include "../../Operators/k_mass/P1ElementwiseKMass.hpp"
#include "../../Operators/k_mass/P1ElementwiseKMassAnnulusMap.hpp"
#include "../../Operators/k_mass/P1ElementwiseKMassIcosahedralShellMap.hpp"
#include "../../Operators/k_p1_coeff_mass/P2ElementwiseKP1CoeffMass.hpp"
#include "../../Operators/k_p1_coeff_mass/P2ElementwiseKP1CoeffMassAnnulusMap.hpp"
#include "../../Operators/k_p1_coeff_mass/P2ElementwiseKP1CoeffMassIcosahedralShellMap.hpp"
#include "../../Operators/mass/P1ElementwiseMass.hpp"
#include "../../Operators/mass/P1ElementwiseMassAnnulusMap.hpp"
#include "../../Operators/mass/P1ElementwiseMassIcosahedralShellMap.hpp"
#include "../../Operators/mass/P2ElementwiseMass.hpp"
#include "../../Operators/mass/P2ElementwiseMassAnnulusMap.hpp"
#include "../../Operators/mass/P2ElementwiseMassIcosahedralShellMap.hpp"
#include "../../Operators/mass_SUPG/P2ElementwiseMassSupg.hpp"
#include "../../Operators/mass_SUPG/P2ElementwiseMassSupgAnnulusMap.hpp"
#include "../../Operators/mass_SUPG/P2ElementwiseMassSupgIcosahedralShellMap.hpp"
#include "../../Operators/shear_heating/P2ElementwiseShearHeating.hpp"
#include "../../Operators/shear_heating/P2ElementwiseShearHeatingAnnulusMap.hpp"
#include "../../Operators/shear_heating/P2ElementwiseShearHeatingIcosahedralShellMap.hpp"
#include "../../Operators/shear_heating_no_surface/P2ElementwiseShearHeatingNoSurface.hpp"
#include "../../Operators/shear_heating_no_surface/P2ElementwiseShearHeatingNoSurfaceAnnulusMap.hpp"
#include "../../Operators/shear_heating_no_surface/P2ElementwiseShearHeatingNoSurfaceIcosahedralShellMap.hpp"
#include "../../Operators/shear_heating_supg/P2ElementwiseShearHeatingSupg.hpp"
#include "../../Operators/shear_heating_supg/P2ElementwiseShearHeatingSupgAnnulusMap.hpp"
#include "../../Operators/shear_heating_supg/P2ElementwiseShearHeatingSupgIcosahedralShellMap.hpp"
#include "../../Operators/shear_heating_supg_no_surface/P2ElementwiseShearHeatingSupgNoSurface.hpp"
#include "../../Operators/shear_heating_supg_no_surface/P2ElementwiseShearHeatingSupgNoSurfaceAnnulusMap.hpp"
#include "../../Operators/shear_heating_supg_no_surface/P2ElementwiseShearHeatingSupgNoSurfaceIcosahedralShellMap.hpp"
#include "../CompositeOperators/P1ToP2GradientOperator.hpp"
#include "../CompositeOperators/P2RhoGMassOperator.hpp"
#include "../CompositeOperators/P2ToP1DivergenceOperator.hpp"
#include "../CompositeOperators/P2ToP1GradRhoRhoDivergenceOperator.hpp"
#include "../CompositeOperators/P2toP1GradRhoRhoOperator.hpp"
#include "../LHS/AdvectionDiffusionOperator.hpp"
#include "../LHS/SaddlePointOperator.hpp"
#include "../RHS/AdvectionDiffusionRHS.hpp"
#include "../RHS/SaddlePointOperatorRHS.hpp"

// Stokes

// #####################
// #### No Operator ####
// #####################

typedef hyteg::NoOperator MC_NoOp;

// #############################
// #### Freeslip Projection ####
// #############################

typedef hyteg::P2ProjectNormalOperator MC_Projection;

// #################
// #### A Block ####
// #################

typedef hyteg::mcoperators::P2VectorElementwiseFullStokesScaledDivdivVectorial                    MC_ABlock;
typedef hyteg::mcoperators::P2VectorElementwiseFullStokesScaledDivdivVectorialAnnulusMap          MC_ABlock_AnnulusMap;
typedef hyteg::mcoperators::P2VectorElementwiseFullStokesScaledDivdivVectorialIcosahedralShellMap MC_ABlock_IcosahedralShellMap;

typedef hyteg::mcoperators::P2VectorElementwiseFullStokesScaledDivdivVectorial           MC_ABlock_Vec;
typedef hyteg::mcoperators::P2VectorElementwiseFullStokesScaledDivdivVectorialAnnulusMap MC_ABlock_Vec_AnnulusMap;
typedef hyteg::mcoperators::P2VectorElementwiseFullStokesScaledDivdivVectorialIcosahedralShellMap
    MC_ABlock_Vec_IcosahedralShellMap;

// ################################
// #### A Block FK Simple Visc ####
// ################################

typedef hyteg::mcoperators::P2VectorElementwiseFullStokesScaledDivdivVectorialFrankKamenetskiiSimpleViscBase
    MC_ABlock_Vec_FK_SimpleVisc;
typedef hyteg::mcoperators::P2VectorElementwiseFullStokesScaledDivdivVectorialFrankKamenetskiiSimpleViscBaseAnnulusMap
    MC_ABlock_Vec_AnnulusMap_FK_SimpleVisc;
typedef hyteg::mcoperators::P2VectorElementwiseFullStokesScaledDivdivVectorialFrankKamenetskiiSimpleViscBaseIcosahedralShellMap
    MC_ABlock_Vec_IcosahedralShellMap_FK_SimpleVisc;

// ##################
// #### BT Block ####
// ##################

typedef hyteg::mcoperators::P1ToP2GradientOperator                    MC_BTBlock;
typedef hyteg::mcoperators::P1ToP2GradientAnnulusMapOperator          MC_BTBlock_AnnulusMap;
typedef hyteg::mcoperators::P1ToP2GradientIcosahedralShellMapOperator MC_BTBlock_IcosahedralShellMap;

// #################
// #### B Block ####
// #################

typedef hyteg::mcoperators::P2ToP1DivergenceOperator                    MC_BBlock;
typedef hyteg::mcoperators::P2ToP1DivergenceAnnulusMapOperator          MC_BBlock_AnnulusMap;
typedef hyteg::mcoperators::P2ToP1DivergenceIcosahedralShellMapOperator MC_BBlock_IcosahedralShellMap;

// #################
// #### ALA B Block ####
// #################

typedef hyteg::mcoperators::P2ToP1DivergenceALAOperator                    MC_BBlock_ALA;
typedef hyteg::mcoperators::P2ToP1DivergenceALAAnnulusMapOperator          MC_BBlock_ALA_AnnulusMap;
typedef hyteg::mcoperators::P2ToP1DivergenceALAIcosahedralShellMapOperator MC_BBlock_ALA_IcosahedralShellMap;

// #######################
// #### Stabilisation ####
// #######################

typedef hyteg::NoOperator MC_Stabilisation;

// #####################################
// #### Temperature to Velocity RHS ####
// #####################################

typedef hyteg::mcoperators::P2RhoGMassOperator                    MC_TemperatureToVelocityRHS;
typedef hyteg::mcoperators::P2RhoGMassAnnulusMapOperator          MC_TemperatureToVelocityRHS_AnnulusMap;
typedef hyteg::mcoperators::P2RhoGMassIcosahedralShellMapOperator MC_TemperatureToVelocityRHS_IcosahedralShellMap;

// ##################################
// #### Velocity to Pressure RHS ####
// ##################################

typedef hyteg::mcoperators::P2toP1GradRhoRhoOperator                    MC_VelocityToPressureRHS;
typedef hyteg::mcoperators::P2toP1GradRhoRhoAnnulusMapOperator          MC_VelocityToPressureRHS_AnnulusMap;
typedef hyteg::mcoperators::P2toP1GradRhoRhoIcosahedralShellMapOperator MC_VelocityToPressureRHS_IcosahedralShellMap;

typedef hyteg::mcoperators::P2toP1GradRhoRhoOperator                    MC_GradRhoRho;
typedef hyteg::mcoperators::P2toP1GradRhoRhoAnnulusMapOperator          MC_GradRhoRho_AnnulusMap;
typedef hyteg::mcoperators::P2toP1GradRhoRhoIcosahedralShellMapOperator MC_GradRhoRho_IcosahedralShellMap;

// Transport

// ##############
// #### Mass ####
// ##############

typedef hyteg::mcoperators::P2ElementwiseMass                    MC_P2Mass;
typedef hyteg::mcoperators::P2ElementwiseMassAnnulusMap          MC_P2Mass_AnnulusMap;
typedef hyteg::mcoperators::P2ElementwiseMassIcosahedralShellMap MC_P2Mass_IcosahedralShellMap;

typedef hyteg::mcoperators::P1ElementwiseMass                    MC_P1Mass;
typedef hyteg::mcoperators::P1ElementwiseMassAnnulusMap          MC_P1Mass_AnnulusMap;
typedef hyteg::mcoperators::P1ElementwiseMassIcosahedralShellMap MC_P1Mass_IcosahedralShellMap;

// ###################
// #### Advection ####
// ###################

typedef hyteg::operatorgeneration::P2ElementwiseAdvection                    MC_Advection;
typedef hyteg::operatorgeneration::P2ElementwiseAdvectionAnnulusMap          MC_Advection_AnnulusMap;
typedef hyteg::operatorgeneration::P2ElementwiseAdvectionIcosahedralShellMap MC_Advection_IcosahedralShellMap;

// ####################
// #### Div K Grad ####
// ####################

typedef hyteg::operatorgeneration::P2ElementwiseDivKGradP1Coefficient                    MC_DivKGrad;
typedef hyteg::operatorgeneration::P2ElementwiseDivKGradP1CoefficientAnnulusMap          MC_DivKGrad_AnnulusMap;
typedef hyteg::operatorgeneration::P2ElementwiseDivKGradP1CoefficientIcosahedralShellMap MC_DivKGrad_IcosahedralShellMap;

typedef hyteg::operatorgeneration::P1ElementwiseDivKGrad                    MC_P1DivKGrad;
typedef hyteg::operatorgeneration::P1ElementwiseDivKGradAnnulusMap          MC_P1DivKGrad_AnnulusMap;
typedef hyteg::operatorgeneration::P1ElementwiseDivKGradIcosahedralShellMap MC_P1DivKGrad_IcosahedralShellMap;

// ##############################
// #### Diffusion Additional ####
// ##############################

typedef hyteg::operatorgeneration::P2ElementwiseDiffusionInvRho                    MC_DiffusionAdditional;
typedef hyteg::operatorgeneration::P2ElementwiseDiffusionInvRhoAnnulusMap          MC_DiffusionAdditional_AnnulusMap;
typedef hyteg::operatorgeneration::P2ElementwiseDiffusionInvRhoIcosahedralShellMap MC_DiffusionAdditional_IcosahedralShellMap;

// ###########################
// #### Adiabatic Heating ####
// ###########################

typedef hyteg::operatorgeneration::P2ElementwiseAdiabaticHeatingAnnulusMap          MC_AdiabaticHeating_AnnulusMap;
typedef hyteg::operatorgeneration::P2ElementwiseAdiabaticHeatingIcosahedralShellMap MC_AdiabaticHeating_IcosahedralShellMap;

// #######################
// #### Shear Heating ####
// #######################

typedef hyteg::mcoperators::P2ElementwiseShearHeating                    MC_ShearHeating;
typedef hyteg::mcoperators::P2ElementwiseShearHeatingAnnulusMap          MC_ShearHeating_AnnulusMap;
typedef hyteg::mcoperators::P2ElementwiseShearHeatingIcosahedralShellMap MC_ShearHeating_IcosahedralShellMap;

// ##################################
// #### Shear Heating No Surface ####
// ##################################

typedef hyteg::mcoperators::P2ElementwiseShearHeatingNoSurface                    MC_ShearHeating_NoSurface;
typedef hyteg::mcoperators::P2ElementwiseShearHeatingNoSurfaceAnnulusMap          MC_ShearHeating_NoSurface_AnnulusMap;
typedef hyteg::mcoperators::P2ElementwiseShearHeatingNoSurfaceIcosahedralShellMap MC_ShearHeating_NoSurface_IcosahedralShellMap;

// ###################
// #### Mass SUPG ####
// ###################

typedef hyteg::mcoperators::P2ElementwiseMassSupg                    MC_P2MassSUPG;
typedef hyteg::mcoperators::P2ElementwiseMassSupgAnnulusMap          MC_P2MassSUPG_AnnulusMap;
typedef hyteg::mcoperators::P2ElementwiseMassSupgIcosahedralShellMap MC_P2MassSUPG_IcosahedralShellMap;

// ########################
// #### Advection SUPG ####
// ########################

typedef hyteg::operatorgeneration::P2ElementwiseSupgAdvection                    MC_AdvectionSUPG;
typedef hyteg::operatorgeneration::P2ElementwiseSupgAdvectionAnnulusMap          MC_AdvectionSUPG_AnnulusMap;
typedef hyteg::operatorgeneration::P2ElementwiseSupgAdvectionIcosahedralShellMap MC_AdvectionSUPG_IcosahedralShellMap;

// ########################
// #### Diffusion SUPG ####
// ########################

typedef hyteg::operatorgeneration::P2ElementwiseSupgDiffusion                    MC_DiffusionSUPG;
typedef hyteg::operatorgeneration::P2ElementwiseSupgDiffusionAnnulusMap          MC_DiffusionSUPG_AnnulusMap;
typedef hyteg::operatorgeneration::P2ElementwiseSupgDiffusionIcosahedralShellMap MC_DiffusionSUPG_IcosahedralShellMap;

// ################################
// #### Adiabatic Heating SUPG ####
// ################################

typedef hyteg::operatorgeneration::P2ElementwiseSupgAdiabaticHeatingAnnulusMap MC_AdiabaticHeatingSUPG_AnnulusMap;
typedef hyteg::operatorgeneration::P2ElementwiseSupgAdiabaticHeatingIcosahedralShellMap
    MC_AdiabaticHeatingSUPG_IcosahedralShellMap;

// ############################
// #### Shear Heating SUPG ####
// ############################

typedef hyteg::mcoperators::P2ElementwiseShearHeatingSupg                    MC_ShearHeatingSUPG;
typedef hyteg::mcoperators::P2ElementwiseShearHeatingSupgAnnulusMap          MC_ShearHeatingSUPG_AnnulusMap;
typedef hyteg::mcoperators::P2ElementwiseShearHeatingSupgIcosahedralShellMap MC_ShearHeatingSUPG_IcosahedralShellMap;

// #######################################
// #### Shear Heating SUPG No Surface ####
// #######################################

typedef hyteg::mcoperators::P2ElementwiseShearHeatingSupgNoSurface           MC_ShearHeatingSUPG_NoSurface;
typedef hyteg::mcoperators::P2ElementwiseShearHeatingSupgNoSurfaceAnnulusMap MC_ShearHeatingSUPG_NoSurface_AnnulusMap;
typedef hyteg::mcoperators::P2ElementwiseShearHeatingSupgNoSurfaceIcosahedralShellMap
    MC_ShearHeatingSUPG_NoSurface_IcosahedralShellMap;

// ####################
// ###### K Mass ######
// ####################

typedef hyteg::mcoperators::P1ElementwiseKMass                    MC_KMass;
typedef hyteg::mcoperators::P1ElementwiseKMassAnnulusMap          MC_KMass_AnnulusMap;
typedef hyteg::mcoperators::P1ElementwiseKMassIcosahedralShellMap MC_KMass_IcosahedralShellMap;

typedef hyteg::mcoperators::P2ElementwiseKP1CoeffMass                    MC_P2KP1Mass;
typedef hyteg::mcoperators::P2ElementwiseKP1CoeffMassAnnulusMap          MC_P2KP1Mass_AnnulusMap;
typedef hyteg::mcoperators::P2ElementwiseKP1CoeffMassIcosahedralShellMap MC_P2KP1Mass_IcosahedralShellMap;

// ####################################
// #### GradRhoRhoDivergence Block ####
// ####################################

typedef hyteg::mcoperators::P2ToP1GradRhoRhoDivergenceOperator                    MC_GradRhoRhoDivergence;
typedef hyteg::mcoperators::P2ToP1GradRhoRhoDivergenceAnnulusMapOperator          MC_GradRhoRhoDivergence_AnnulusMap;
typedef hyteg::mcoperators::P2ToP1GradRhoRhoDivergenceIcosahedralShellMapOperator MC_GradRhoRhoDivergence_IcosahedralShellMap;