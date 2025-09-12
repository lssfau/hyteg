/*
 * Copyright (c) 2025 Andreas Burkhart.
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
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionHistory.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2InjectionRestriction.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/numerictools/BDFScheme.hpp"
#include "hyteg/numerictools/CrankNicolsonScheme.hpp"
#include "hyteg/operators/GEMV.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/solvers/preconditioners/IdentityPreconditioner.hpp"

#include "../../Utility/Data/DataLoader.hpp"
#include "../../Utility/Density/ConstantDensity.hpp"
#include "../../Utility/Density/ExponentialDensity.hpp"
#include "../../Utility/Density/LinearDensityInterpolation.hpp"
#include "../../Utility/Density/PressureProfileDensityBilinearInterpolation.hpp"
#include "../../Utility/LHS/AdvectionDiffusionOperator.hpp"
#include "../../Utility/LHS/SaddlePointOperator.hpp"
#include "../../Utility/OperatorTools/CombinedABlockOperator.hpp"
#include "../../Utility/OperatorTools/OptimisationTypedefs.hpp"
#include "../../Utility/Parameters/NondimensionalisationParameters.hpp"
#include "../../Utility/Pressure/ConstantPressure.hpp"
#include "../../Utility/Pressure/LinearPressureInterpolation.hpp"
#include "../../Utility/RHS/SaddlePointOperatorRHS.hpp"
#include "../../Utility/Solver/ABlock/ABlockCGCoarseGridSolver.hpp"
#include "../../Utility/Solver/ABlock/ABlockCGOuterLoopSolver.hpp"
#include "../../Utility/Solver/ABlock/ABlockCGVectorMassSolver.hpp"
#include "../../Utility/Solver/ABlock/ABlockChebyshevSmoother.hpp"
#include "../../Utility/Solver/ABlock/ABlockInverseDiagSolver.hpp"
#include "../../Utility/Solver/ABlock/ABlockMinResCoarseGridSolver.hpp"
#include "../../Utility/Solver/ABlock/ABlockMinResOuterLoopSolver.hpp"
#include "../../Utility/Solver/ABlock/ABlockMultigridSolver.hpp"
#include "../../Utility/Solver/ABlock/ABlockPETScCoarseGridSolver.hpp"
#include "../../Utility/Solver/AdvectionDiffusion/AdvectionDiffusionCGSPDPreconditioner.hpp"
#include "../../Utility/Solver/AdvectionDiffusion/AdvectionDiffusionFGMRESSolver.hpp"
#include "../../Utility/Solver/SaddlePoint/SaddlePointAdjointInexactUzawaSmoother.hpp"
#include "../../Utility/Solver/SaddlePoint/SaddlePointBlockApproximationFactorisationSmoother.hpp"
#include "../../Utility/Solver/SaddlePoint/SaddlePointFGMRESSolver.hpp"
#include "../../Utility/Solver/SaddlePoint/SaddlePointInexactUzawaSmoother.hpp"
#include "../../Utility/Solver/SaddlePoint/SaddlePointSymmetricUzawaSmoother.hpp"
#include "../../Utility/Solver/Schur/SchurBFBTSolver.hpp"
#include "../../Utility/Solver/Schur/SchurCGCoarseGridSolver.hpp"
#include "../../Utility/Solver/Schur/SchurCGMassSolver.hpp"
#include "../../Utility/Solver/Schur/SchurCGOuterLoopSolver.hpp"
#include "../../Utility/Solver/Schur/SchurCGPoissonNeumannSolver.hpp"
#include "../../Utility/Solver/Schur/SchurChebyshevSmoother.hpp"
#include "../../Utility/Solver/Schur/SchurInverseDiagMassSolver.hpp"
#include "../../Utility/Solver/Schur/SchurInverseLumpedMassSolver.hpp"
#include "../../Utility/Solver/Schur/SchurMultigridSolver.hpp"
#include "../../Utility/Solver/Schur/SchurPoissonNeumannChebyshevSmoother.hpp"
#include "../../Utility/Solver/Schur/SchurPoissonNeumannMultigridSolver.hpp"
#include "../../Utility/Temperature/ConstantTemperature.hpp"
#include "../../Utility/Temperature/ExponentialTemperature.hpp"
#include "../../Utility/Temperature/LinearTemperature.hpp"
#include "../../Utility/Temperature/LinearTemperatureInterpolation.hpp"
#include "../../Utility/Temperature/RelativeRandomTemperature.hpp"
#include "../../Utility/Temperature/SphericalHarmonicsTemperature.hpp"
#include "../../Utility/Viscosity/ArrheniusViscosity.hpp"
#include "../../Utility/Viscosity/ConstantViscosity.hpp"
#include "../../Utility/Viscosity/ExponentialViscosity.hpp"
#include "../../Utility/Viscosity/LinearViscosityInterpolation.hpp"
#include "../../Utility/Viscosity/SimpleSpaceDependentViscosityProfileWithJumps.hpp"
#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "elementwise_dof_value_operator/generated/p2_shear_heat_T_p2_dep_eta_blending_q6_ElementwiseOperator.hpp"
#include "terraneo/helpers/ConvectionToolbox.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;
using namespace hyteg::convectionToolbox;
using namespace MantleConvection;

// #define NoBlend

#ifdef NoBlend
// Define operator types
typedef MantleConvection::AdvectionDiffusionOperator< MC_P2Mass,
                                                      MC_Advection,
                                                      MC_NoOp,
                                                      MC_NoOp,
                                                      MC_NoOp,

                                                      MC_NoOp,
                                                      MC_NoOp,
                                                      MC_NoOp,
                                                      MC_NoOp >
    transportLHSType;

typedef MantleConvection::AdvectionDiffusionOperatorRHS< MC_P2Mass,
                                                         MC_NoOp,
                                                         MC_NoOp,
                                                         MC_NoOp,

                                                         MC_NoOp,
                                                         MC_NoOp,
                                                         MC_NoOp,
                                                         MC_NoOp >
    transportRHSType;
#else
// Define operator types
typedef MantleConvection::AdvectionDiffusionOperator< MC_P2Mass_AnnulusMap,
                                                      MC_Advection_AnnulusMap,
                                                      MC_NoOp,
                                                      MC_NoOp,
                                                      MC_NoOp,

                                                      MC_NoOp,
                                                      MC_NoOp,
                                                      MC_NoOp,
                                                      MC_NoOp >
    transportLHSType;

typedef MantleConvection::AdvectionDiffusionOperatorRHS< MC_P2Mass_AnnulusMap,
                                                         MC_NoOp,
                                                         MC_NoOp,
                                                         MC_NoOp,

                                                         MC_NoOp,
                                                         MC_NoOp,
                                                         MC_NoOp,
                                                         MC_NoOp >
    transportRHSType;
#endif

int main( int argc, char** argv )
{
   // Create walberla & MPI environment
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   // PETSc
   PETScManager petscManager( &argc, &argv );
#endif

   std::shared_ptr< walberla::config::Config > cfg_ = std::make_shared< walberla::config::Config >();
   cfg_->readParameterFile( "./MC_PureAdvectionConvTest.prm" );
   walberla::config::Config::BlockHandle parameters_ = cfg_->getOneBlock( "Parameters" );

   // parameters
   real_t delta_t   = parameters_.getParameter< real_t >( "delta_t" );
   bool   vtkOutput = parameters_.getParameter< bool >( "vtk" );
   ;
   real_t t_init    = parameters_.getParameter< real_t >( "t_init" );
   real_t t_current = t_init;
   real_t run_time  = parameters_.getParameter< real_t >( "run_time" );
   uint_t steps     = (uint_t) std::ceil( run_time / delta_t );
   bool   MMOC      = parameters_.getParameter< bool >( "MMOC" );
#ifdef NoBlend
   bool blending = false;
#else
   bool blending          = true;
#endif
   real_t rel_tolerance = parameters_.getParameter< real_t >( "rel_tolerance" );
   real_t abs_tolerance = parameters_.getParameter< real_t >( "abs_tolerance" );

   uint_t minLevel_      = parameters_.getParameter< uint_t >( "minLevel" );
   uint_t maxLevel_      = parameters_.getParameter< uint_t >( "maxLevel" );
   int    BDFOrder       = parameters_.getParameter< int >( "BDFOrder" );
   bool   lowMemoryMode_ = false;

//storage
#ifdef NoBlend
   uint_t n_       = 8;
   auto   storage_ = createRectangleStorage( n_, n_, Point2D( { 0.5, 0.5 } ), Point2D( { 1.5, 1.5 } ) );
   auto   BC_      = createRectangleBoundaryConditions(
       DirichletBoundary, DirichletBoundary, DirichletBoundary, DirichletBoundary, DirichletBoundary );
#else
   auto storage_          = createAnnulusStorage( 24, 8, 0.5, 1.5, blending );
   auto BC_               = createAnnulusBoundaryConditions( DirichletBoundary, DirichletBoundary );
#endif

   // constants
   const real_t k_     = parameters_.getParameter< real_t >( "k" );
   const real_t C_p_   = 1.0;
   const real_t Pe_    = 1.0;
   const real_t H_     = 1.0;
   const real_t Ra_    = 1.0;
   const real_t Di_    = 1.0;
   const real_t alpha_ = 1.0;

   // define std::functions
   std::function< real_t( const hyteg::Point3D&, real_t ) > eta_temp = [&]( const hyteg::Point3D& x, real_t temp ) {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( temp );
      return 1;
   };

   std::function< real_t( const hyteg::Point3D&, real_t ) > eta = [&]( const hyteg::Point3D& x, real_t t ) {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( t );
      return 1;
   };

   std::function< real_t( const hyteg::Point3D& ) > rho = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return 1;
   };

   std::function< real_t( const hyteg::Point3D&, real_t ) > T_fct = [&]( const hyteg::Point3D& x, real_t t ) {
      WALBERLA_UNUSED( x );
      return exp( ( 1.0 / 4.0 ) *
                  ( -( ( x[0] + sin( t ) ) * ( x[0] + sin( t ) ) ) - ( ( x[1] - cos( t ) ) * ( x[1] - cos( t ) ) ) ) / k_ );
   };
   std::function< real_t( const hyteg::Point3D&, real_t ) > uX = [&]( const hyteg::Point3D& x, real_t t ) {
      WALBERLA_UNUSED( x );
      return -x[1];
   };
   std::function< real_t( const hyteg::Point3D&, real_t ) > uY = [&]( const hyteg::Point3D& x, real_t t ) {
      WALBERLA_UNUSED( x );
      return x[0];
   };
   std::function< real_t( const hyteg::Point3D&, real_t ) > uZ = [&]( const hyteg::Point3D& x, real_t t ) {
      WALBERLA_UNUSED( x );
      return 0;
   };

   std::function< real_t( const hyteg::Point3D&, real_t ) > ctrl_RHS = [&]( const hyteg::Point3D& x, real_t t ) {
      WALBERLA_UNUSED( x );
      return 0;
   };

   std::function< real_t( const hyteg::Point3D& ) > gX = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return -x[0];
   };
   std::function< real_t( const hyteg::Point3D& ) > gY = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return -x[1];
   };
   std::function< real_t( const hyteg::Point3D& ) > gZ = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return 0;
   };

   // ############################
   // ### Define FEM functions ###
   // ############################
   uint_t numberOfStatesT = (uint_t) BDFOrder + 2;
   uint_t numberOfStatesU = (uint_t) BDFOrder + 2;

   auto T_ = std::make_shared< FunctionHistory< P2Function< real_t >, real_t > >( numberOfStatesT, minLevel_, maxLevel_ );
   for ( uint_t i = 0; i < T_->getMemoryCapacity(); i++ )
   {
      std::stringstream fName;
      fName << "T history " << i;
      auto TFct = createSharedP2TemperatureFunction( fName.str(), storage_, minLevel_, maxLevel_, BC_ );
      T_->addFunction( TFct );
   }

   auto u_ = std::make_shared< FunctionHistory< P2VectorFunction< real_t >, real_t > >( numberOfStatesU, minLevel_, maxLevel_ );
   for ( uint_t i = 0; i < u_->getMemoryCapacity(); i++ )
   {
      std::stringstream fName;
      fName << "u history " << i;
      auto upFct = createSharedP2VectorVelocityFunction( fName.str(), storage_, minLevel_, maxLevel_, BC_ );
      u_->addFunction( upFct );
   }

   auto T_extra_ = createSharedP2TemperatureFunction( "T_extra", storage_, minLevel_, maxLevel_, BC_ );
   auto u_extra_ = createSharedP2VectorVelocityFunction( "u_extra", storage_, minLevel_, maxLevel_, BC_ );
   auto T_ref    = createSharedP2TemperatureFunction( "T_ref", storage_, minLevel_, maxLevel_, BC_ );
   auto tmp      = createSharedP2TemperatureFunction( "tmp", storage_, minLevel_, maxLevel_, BC_ );

   auto rho_ =
       std::make_shared< P1Function< real_t > >( "rho", storage_, minLevel_, maxLevel_, BoundaryCondition::createAllInnerBC() );
   auto inv_rho_ = std::make_shared< P1Function< real_t > >(
       "inv_rho", storage_, minLevel_, maxLevel_, BoundaryCondition::createAllInnerBC() );

   auto rhs_     = createSharedP2TemperatureFunction( "rhs", storage_, minLevel_, maxLevel_, BC_ );
   auto rhsLast_ = createSharedP2TemperatureFunction( "rhsLast", storage_, minLevel_, maxLevel_, BC_ );

   auto eta_ref_ = std::make_shared< P1Function< real_t > >(
       "eta ref", storage_, minLevel_, maxLevel_, BoundaryCondition::createAllInnerBC() );
   auto eta_extra_ = std::make_shared< P1Function< real_t > >(
       "eta extra", storage_, minLevel_, maxLevel_, BoundaryCondition::createAllInnerBC() );

   // ###########################################
   // ### Define advection diffusion operator ###
   // ###########################################

   // auto timeScheme =
   //     std::make_shared< CrankNicolsonScheme< P2Function< real_t >, transportLHSType::AdditiveMassType > >( rhsLast_ );

   auto timeScheme = std::make_shared< BDFScheme< P2Function< real_t >, transportLHSType::AdditiveMassType > >( BDFOrder );

   auto MassOperator = std::make_shared< transportLHSType::MassOperatorTypeInternal >( storage_, minLevel_, maxLevel_ );

#ifdef NoBlend
   auto AdvectionOperator = std::make_shared< transportLHSType::AdvectionOperatorTypeInternal >(
       storage_, minLevel_, maxLevel_, u_extra_->component( 0 ), u_extra_->component( 1 ), u_extra_->component( 0 ) );
#else
   auto AdvectionOperator = std::make_shared< transportLHSType::AdvectionOperatorTypeInternal >(
       storage_, minLevel_, maxLevel_, u_extra_->component( 0 ), u_extra_->component( 1 ) );
#endif

   real_t AdvectionScaling           = MMOC ? real_c( 0 ) : real_c( 1 );
   real_t DiffusionScaling           = real_c( 0 );
   real_t DiffusionAdditionalScaling = real_c( 0 );
   real_t AdiabaticHeatingScaling    = real_c( 0 );

   auto transportOperator_ = std::make_shared< transportLHSType >( storage_,
                                                                   minLevel_,
                                                                   maxLevel_,
                                                                   T_,
                                                                   MassOperator,
                                                                   AdvectionOperator,
                                                                   nullptr,
                                                                   nullptr,
                                                                   nullptr,
                                                                   AdvectionScaling,
                                                                   DiffusionScaling,
                                                                   DiffusionAdditionalScaling,
                                                                   AdiabaticHeatingScaling,
                                                                   timeScheme );

   // ###############################################
   // ### Define advection diffusion RHS operator ###
   // ###############################################
   std::function< real_t( const hyteg::Point3D&, real_t ) > eta_extra_current = [&]( const hyteg::Point3D& x, real_t temp ) {
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( temp );
      return eta( x, t_current );
   };

   std::function< real_t( const Point3D& ) > ShearHeatScalingFct = [&]( const hyteg::Point3D& x ) {
      return real_c( 1 ) / rho( x );
   };

   auto ShearHeatingOperator =
       std::make_shared< p2_shear_heat_T_p2_dep_eta_blending_q6_ElementwiseOperator >( storage_,
                                                                                       minLevel_,
                                                                                       maxLevel_,
                                                                                       u_extra_->component( 0 ),
                                                                                       u_extra_->component( 1 ),
                                                                                       u_extra_->component( 0 ),
                                                                                       *T_extra_,
                                                                                       eta_temp,
                                                                                       ShearHeatScalingFct );

   // auto ShearHeatingOperator = std::make_shared< transportRHSType::ShearHeatingOperatorTypeInternal >(
   //     storage_, minLevel_, maxLevel_, *eta_extra_, *inv_rho_, u_extra_->component( 0 ), u_extra_->component( 1 ) );

   // auto ShearHeatingOperatorCN = std::make_shared< transportRHSType::ShearHeatingOperatorTypeInternal >(
   //     storage_, minLevel_, maxLevel_, *eta_extra_, *inv_rho_, u_->getState(1).component(0), u_->getState(1).component( 1 ) );

   real_t InternalHeatingScaling = real_c( 0 );
   real_t ShearHeatingScaling    = real_c( 0 );

   // InternalHeatingScaling = 0;
   // ShearHeatingScaling = 0;

   auto transportOperator_RHS_ = std::make_shared< transportRHSType >( storage_,
                                                                       minLevel_,
                                                                       maxLevel_,
                                                                       T_,
                                                                       MassOperator,
                                                                       nullptr,
                                                                       nullptr,
                                                                       nullptr,
                                                                       InternalHeatingScaling,
                                                                       ShearHeatingScaling,
                                                                       AdiabaticHeatingScaling,
                                                                       timeScheme,
                                                                       nullptr,
                                                                       nullptr,
                                                                       nullptr,
                                                                       nullptr,
                                                                       real_c( 0 ),
                                                                       real_c( 0 ),
                                                                       real_c( 0 ) );

   // ###########################
   // ### Test Implementation ###
   // ###########################

   P2ElementwiseBlendingMassOperator P2Mass( storage_, minLevel_, maxLevel_ );

   // init MMOC
   auto mmocTransport =
       std::make_shared< MMOCTransport< P2Function< real_t > > >( storage_, minLevel_, maxLevel_, TimeSteppingScheme::RK4 );

   // interpolate
   std::function< real_t( const hyteg::Point3D&, const std::vector< real_t >& ) > eta_extra_temp =
       [&]( const hyteg::Point3D& x, const std::vector< real_t >& temp ) {
          WALBERLA_UNUSED( x );
          WALBERLA_UNUSED( temp );
          return eta_temp( x, temp[0] );
       };

   std::function< real_t( const hyteg::Point3D& ) > T_current = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return T_fct( x, t_current );
   };
   std::function< real_t( const hyteg::Point3D& ) > ctrl_RHS_current = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return ctrl_RHS( x, t_current );
   };
   // std::function< real_t( const hyteg::Point3D& ) > ctrl_RHS_MMOC_current = [&]( const hyteg::Point3D& x ) {
   //    WALBERLA_UNUSED( x );
   //    return ctrl_RHS_MMOC( x, t_current );
   // };
   std::function< real_t( const hyteg::Point3D& ) > uX_current = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return uX( x, t_current );
   };
   std::function< real_t( const hyteg::Point3D& ) > uY_current = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return uY( x, t_current );
   };
   std::function< real_t( const hyteg::Point3D& ) > uZ_current = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return uZ( x, t_current );
   };
   std::function< real_t( const hyteg::Point3D& ) > eta_current = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return eta( x, t_current );
   };

   std::function< real_t( const hyteg::Point3D& ) > T_Minus1 = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return T_fct( x, t_current - delta_t );
   };
   std::function< real_t( const hyteg::Point3D& ) > T_0 = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return T_fct( x, t_init );
   };
   std::function< real_t( const hyteg::Point3D& ) > uX_Minus1 = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return uX( x, t_current - delta_t );
   };
   std::function< real_t( const hyteg::Point3D& ) > uX_0 = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return uX( x, t_init );
   };
   std::function< real_t( const hyteg::Point3D& ) > uY_Minus1 = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return uY( x, t_current - delta_t );
   };
   std::function< real_t( const hyteg::Point3D& ) > uY_0 = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return uY( x, t_init );
   };
   std::function< real_t( const hyteg::Point3D& ) > uZ_Minus1 = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return uZ( x, t_current - delta_t );
   };
   std::function< real_t( const hyteg::Point3D& ) > uZ_0 = [&]( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return uZ( x, t_init );
   };

   // past states
   t_current = t_init - BDFOrder * delta_t;

   for ( uint_t i = 0; i < BDFOrder; i++ )
   {
      T_->newState( real_c( delta_t ) );
      u_->newState( real_c( delta_t ) );
      T_->getState( 0 ).interpolate( T_current, maxLevel_, All );
      u_->getState( 0 ).interpolate( { uX_current, uY_current, uZ_current }, maxLevel_, All );

      t_current += delta_t;
   }

   t_current = t_init;

   T_->newState( real_c( delta_t ) );
   u_->newState( real_c( delta_t ) );
   T_->getState( 0 ).interpolate( T_current, maxLevel_, All );
   u_->getState( 0 ).interpolate( { uX_current, uY_current, uZ_current }, maxLevel_, All );

   rho_->interpolate( rho, maxLevel_, All );
   inv_rho_->assign( { real_c( 1 ) }, { *rho_ }, maxLevel_, All );
   inv_rho_->invertElementwise( maxLevel_, hyteg::All, false );
   eta_ref_->interpolate( eta_current, maxLevel_, All );
   eta_extra_->interpolate( eta_current, maxLevel_, All );

   {
      tmp->interpolate( ctrl_RHS_current, maxLevel_, All );
      u_extra_->assign( { real_c( 1 ) }, { u_->getState( 0 ) }, maxLevel_, All );
      T_extra_->assign( { real_c( 1 ) }, { T_->getState( 0 ) }, maxLevel_, All );
      eta_extra_->interpolate( eta_current, maxLevel_, All );
      transportOperator_->applyTimeIndependent( T_->getState( 0 ), *rhsLast_, maxLevel_, All, Replace );
      rhsLast_->assign( { real_c( -1 ) }, { *rhsLast_ }, maxLevel_, All );
      P2Mass.apply( *tmp, *rhsLast_, maxLevel_, All, Add );
      transportOperator_RHS_->applyTimeIndependent( T_->getState( 0 ), *rhsLast_, maxLevel_, All, Add );
   }

   // iteration callback
   // real_t initRes = 1.0;
   // real_t initErr = 1.0;

   auto stopIterationCallback = [&]( uint_t                      iteration,
                                     const transportLHSType&     _A,
                                     const P2Function< real_t >& _u,
                                     const P2Function< real_t >& _b,
                                     uint_t                      _level,
                                     real_t                      ApproxError,
                                     real_t                      wNorm ) {
      WALBERLA_UNUSED( _A );
      WALBERLA_UNUSED( _u );
      WALBERLA_UNUSED( _b );

      // // check residual
      // _A.apply( _u, *tmp, _level, hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary );
      // tmp->assign( { 1.0, -1.0 }, { *tmp, _b }, _level, hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary );
      // real_t tempResidual =
      //     std::sqrt( tmp->dotGlobal( *tmp, _level, hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary ) );

      // if ( iteration == 0 )
      // {
      //    initRes = tempResidual;
      // }
      // real_t res_rel = tempResidual / initRes;

      // //check error
      // tmp->assign( { 1.0, -1.0 }, { _u, *T_ref }, _level, All );

      // P2Function< real_t > Merr( "Merr", storage_, _level, _level );

      // P2Mass.apply( *tmp, Merr, _level, All );

      // real_t discr_l2_err =
      //     std::sqrt( tmp->dotGlobal( Merr, _level, hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary ) );

      // if ( iteration == 0 )
      // {
      //    initErr = discr_l2_err;
      // }

      // real_t discr_l2_err_rel = discr_l2_err / initErr;

      // WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
      //     "[Iteration Callback] iter %3d | abs. res: %10.5e | rel. res: %10.5e | L2 Error: %10.5e | L2 Error rel: %10.5e",
      //     iteration,
      //     tempResidual,
      //     res_rel,
      //     discr_l2_err,
      //     discr_l2_err_rel ) );

      // if ( res_rel < rel_tolerance )
      // {
      //    WALBERLA_LOG_INFO_ON_ROOT( "relative residual threshold of " << std::scientific << rel_tolerance
      //                                                                 << " reached. Finishing solver loop." );
      //    return true;
      // }

      return false;
   };

   // solver
   // Advection Diffusion solver
   auto AdvectionDiffusionSPDSubOperator_ = std::make_shared< AdvectionDiffusionSPDSubOperator< transportLHSType > >(
       storage_, minLevel_, maxLevel_, transportOperator_ );

   auto AdvectionDiffusionCGSPDPreconditioner_ = std::make_shared<
       hyteg::CGSolver< AdvectionDiffusionSPDSubOperator< transportLHSType > > >(
       storage_,
       minLevel_,
       maxLevel_,
       1000,
       rel_tolerance,
       abs_tolerance,
       std::make_shared< AdvectionDiffusionIdentityPreconditioner< AdvectionDiffusionSPDSubOperator< transportLHSType > > >() );
   AdvectionDiffusionCGSPDPreconditioner_->setPrintInfo( false );

   auto tempSubstituteSolver_ =
       std::make_shared< SubstitutePreconditioner< transportLHSType, AdvectionDiffusionSPDSubOperator< transportLHSType > > >(
           AdvectionDiffusionCGSPDPreconditioner_, AdvectionDiffusionSPDSubOperator_ );

   auto tempSolver_ = std::make_shared< FGMRESSolver< transportLHSType > >( storage_,
                                                                            minLevel_,
                                                                            maxLevel_,
                                                                            1000,
                                                                            rel_tolerance,
                                                                            abs_tolerance,
                                                                            tempSubstituteSolver_,
                                                                            1000,
                                                                            1e-16,
                                                                            0.0,
                                                                            lowMemoryMode_,
                                                                            stopIterationCallback,
                                                                            false );
   tempSolver_->setPrintInfo( true );

   for ( uint_t iter = 0; iter < steps; iter++ )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "step " << iter << "/" << steps );

      // extrapolate
      // u_->extrapolate( 1, delta_t, *u_extra_, maxLevel_, hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary );
      // T_->extrapolate( 1, delta_t, *T_extra_, maxLevel_, hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary );

      // u_->extrapolate( 1, delta_t, *u_extra_, maxLevel_, hyteg::All );
      // T_->extrapolate( 1, delta_t, *T_extra_, maxLevel_, hyteg::All );

      // eta_extra_->interpolate( eta_extra, { T_extra_->getVertexDoFFunction() }, maxLevel_, All );

      t_current += delta_t;

      // T_extra_->interpolate( T_current, maxLevel_, hyteg::DirichletBoundary );
      // u_extra_->interpolate( { uX_current, uY_current, uZ_current }, maxLevel_, hyteg::DirichletBoundary );

      T_extra_->interpolate( T_current, maxLevel_, hyteg::All );
      u_extra_->interpolate( { uX_current, uY_current, uZ_current }, maxLevel_, hyteg::All );
      eta_extra_->interpolate( eta_current, maxLevel_, All );

      T_->newState( delta_t );
      u_->newState( delta_t );

      u_->getState( 0 ).interpolate( { uX_current, uY_current, uZ_current }, maxLevel_, All );
      T_ref->interpolate( T_current, maxLevel_, All );
      eta_ref_->interpolate( eta_current, maxLevel_, All );

      // mmoc transport
      std::shared_ptr< P2Function< real_t > > tmpTempStep;
      if ( MMOC )
      {
         tmpTempStep = getTemporaryFunction< P2Function< real_t > >( storage_, minLevel_, maxLevel_ );
         tmpTempStep->setBoundaryCondition( BC_.bcTemperature_ );

         tmpTempStep->assign( { real_c( 1 ) }, { T_->getState( 1 ) }, maxLevel_, All );

         WALBERLA_LOG_INFO_ON_ROOT( "Perform MMOC advection step..." );
         mmocTransport->step( T_->getState( 1 ), *u_extra_, *u_extra_, maxLevel_, hyteg::All, delta_t, 1, true );

         T_->getState( 1 ).interpolate( T_current, maxLevel_, DirichletBoundary );

         T_->getState( 0 ).assign( { real_c( 1 ) }, { T_->getState( 1 ) }, maxLevel_, All );
      }
      else
      {
         tmp->interpolate( ctrl_RHS_current, maxLevel_, All );
         P2Mass.apply( *tmp, *rhs_, maxLevel_, All );
         transportOperator_RHS_->apply( *T_extra_, *rhs_, maxLevel_, All, Add );
      }

      if ( MMOC )
      {
         // restore old state 1 in case of mmoc
         T_->getState( 1 ).assign( { real_c( 1 ) }, { *tmpTempStep }, maxLevel_, All );
      }
      else
      {
         T_->getState( 0 ).interpolate( T_current, maxLevel_, DirichletBoundary );
         stopIterationCallback( 0, *transportOperator_, T_->getState( 0 ), *rhs_, maxLevel_, 0, 0 );
         tempSolver_->solve( *transportOperator_, T_->getState( 0 ), *rhs_, maxLevel_ );
      }

      {
         tmp->interpolate( ctrl_RHS_current, maxLevel_, All );
         u_extra_->assign( { real_c( 1 ) }, { u_->getState( 0 ) }, maxLevel_, All );
         T_extra_->assign( { real_c( 1 ) }, { T_->getState( 0 ) }, maxLevel_, All );
         eta_extra_->interpolate( eta_extra_temp, { T_extra_->getVertexDoFFunction() }, maxLevel_, All );
         transportOperator_->applyTimeIndependent( T_->getState( 0 ), *rhsLast_, maxLevel_, All, Replace );
         rhsLast_->assign( { real_c( -1 ) }, { *rhsLast_ }, maxLevel_, All );
         P2Mass.apply( *tmp, *rhsLast_, maxLevel_, All, Add );
         transportOperator_RHS_->applyTimeIndependent( T_->getState( 0 ), *rhsLast_, maxLevel_, All, Add );
      }

      // output
      if ( vtkOutput )
      {
         auto vtkOutput_ = std::make_shared< VTKOutput >( ".", "MC_PureAdvectionConvTest", storage_, 1 );

         for ( uint_t dimension = 0; dimension < ( storage_->hasGlobalCells() ? 3 : 2 ); dimension++ )
         {
            vtkOutput_->add( u_->getFunctionByIndex( 0 )[dimension] );
         }

         vtkOutput_->add( T_->getFunctionByIndex( 0 ) );
         vtkOutput_->add( *rhs_ );

         vtkOutput_->add( *T_extra_ );
         vtkOutput_->add( *u_extra_ );
         vtkOutput_->add( *tmp );

         vtkOutput_->add( *eta_extra_ );
         vtkOutput_->add( *eta_ref_ );

         vtkOutput_->add( *rho_ );
         vtkOutput_->add( *T_ref );

         vtkOutput_->write( maxLevel_, iter );
      }
   }

   tmp->assign( { 1.0, -1.0 }, { T_->getState( 0 ), *T_ref }, maxLevel_, All );

   P2Function< real_t > Merr( "Merr", storage_, maxLevel_, maxLevel_ );

   P2Mass.apply( *tmp, Merr, maxLevel_, All );

   real_t discr_l2_err =
       std::sqrt( tmp->dotGlobal( Merr, maxLevel_, hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary ) );
   WALBERLA_LOG_INFO_ON_ROOT( "finalAbsL2Error = " << discr_l2_err );
}