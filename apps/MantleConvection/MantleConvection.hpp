/*
 * Copyright (c) 2023-2025 Andreas Burkhart.
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

#include <cmath>
#include <filesystem>
#include <fstream>

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2InjectionRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/indexing/ConsistentEnumeration.hpp"
#include "hyteg/memory/TempFunctionManager.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/numerictools/CFDHelpers.hpp"
#include "hyteg/petsc/PETScHDF5FunctionSave.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/preconditioners/IdentityPreconditioner.hpp"
#include "hyteg/sparseassembly/FileWritingVector.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

#include "Utility/Density/TemperatureDependentDensityModel.hpp"
#include "Utility/OperatorTools/OperatorUpdater.hpp"
#include "Utility/OperatorTools/OperatorTypedefs.hpp"
#include "Utility/Parameters/NondimensionalisationParameters.hpp"
#include "Utility/Pressure/PressureModel.hpp"
#include "Utility/Solver/ABlock/ABlockSolver.hpp"
#include "Utility/Solver/AdvectionDiffusion/AdvectionDiffusionSolver.hpp"
#include "Utility/Solver/SaddlePoint/SaddlePointSolver.hpp"
#include "Utility/Solver/Schur/SchurSolver.hpp"
#include "Utility/Temperature/TemperatureModel.hpp"
#include "Utility/Viscosity/TemperatureDependentViscosityModel.hpp"
#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "terraneo/helpers/ConvectionToolbox.hpp"
#include "terraneo/plates/PlateVelocityProvider.hpp"
#include "terraneo/sphericalharmonics/SphericalHarmonicsTool.hpp"

#ifdef HYTEG_BUILD_WITH_ADIOS2
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointExporter.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointImporter.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#endif

using terraneo::SphericalHarmonicsTool;
using walberla::math::pi;
using namespace hyteg;
using namespace hyteg::convectionToolbox;

namespace MantleConvection {

template < class SaddlePointOperatorType_,
           class SaddlePointRHSOperatorType_,
           class AdvectionDiffusionOperatorType_,
           class AdvectionDiffusionRHSOperatorType_,
           class TemperatureMassOperatorType_        = MC_P1Mass_IcosahedralShellMap,
           class TemperatureRestrictionOperatorType_ = hyteg::P1toP1LinearRestriction< real_t >,
           class DensityFunctionType_                = hyteg::P1Function< real_t >,
           class ViscosityFunctionType_              = hyteg::P1Function< real_t >,
           bool CalculateInvRhoAndEta_               = true,
           bool InterpolateDensity_                  = true,
           bool InterpolateViscosity_                = true >
class MantleConvectionModel
{
 public:
   typedef SaddlePointOperatorType_            SaddlePointOperatorType;
   typedef SaddlePointRHSOperatorType_         SaddlePointRHSOperatorType;
   typedef AdvectionDiffusionOperatorType_     AdvectionDiffusionOperatorType;
   typedef AdvectionDiffusionRHSOperatorType_  AdvectionDiffusionRHSOperatorType;
   typedef TemperatureMassOperatorType_        TemperatureMassOperatorType;
   typedef TemperatureRestrictionOperatorType_ TemperatureRestrictionOperatorType;
   typedef DensityFunctionType_                DensityFunctionType;
   typedef ViscosityFunctionType_              ViscosityFunctionType;

   MantleConvectionModel( uint_t dim, std::string parameterfile = "./parameters_MC.prm", std::string prefix = "" )
   : MantleConvectionModel( dim, parameterfile, "", prefix )
   {}

   MantleConvectionModel( uint_t dim, std::string parameterfile, std::string setupPrimitiveStorageFile, std::string prefix )
   : rank_( uint_c( walberla::mpi::MPIManager::instance()->rank() ) )
   , size_( uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) )
   , prefix_( prefix )
   , blending_( true )
   , saveNextCheckpoint_( false )
   , saveNextOutput_( false )
   , dim_( dim )
   , solverFlag_( hyteg::Inner | hyteg::NeumannBoundary | hyteg::FreeslipBoundary )
   , dirichletFlag_( hyteg::DirichletBoundary )
   , allFlag_( hyteg::All )
   , velocityProjectionFlag_( hyteg::FreeslipBoundary )
   , stepCounter_( 0 )
   , dt_( real_c( 0 ) )
   , currentTime_( real_c( 0 ) )
   , secondsPerMillionYears_( real_c( 3.1536e13 ) )
   , currentPlateTime_( real_c( 0 ) )
   , loadStorageFile_( !setupPrimitiveStorageFile.empty() )
   , fileSetupPrimitiveStorage_( setupPrimitiveStorageFile )
   {
      // load parameters
      cfg_ = std::make_shared< walberla::config::Config >();
      cfg_->readParameterFile( parameterfile.c_str() );
      parameters_ = cfg_->getOneBlock( "Parameters" );

      // clang-format off
      minLevel_                                            = parameters_.getParameter< uint_t >( prefix_ + std::string( "minLevel" ) );
      maxLevel_                                            = parameters_.getParameter< uint_t >( prefix_ + std::string( "maxLevel" ) );
      writeFrequencyOutput_                                = parameters_.getParameter< uint_t >( prefix_ + std::string( "writeFrequencyOutput" ) );
      writeFrequencyCheckpoint_                            = parameters_.getParameter< uint_t >( prefix_ + std::string( "writeFrequencyCheckpoint" ) );
      surfaceBoundaryType_                                 = parameters_.getParameter< uint_t >( prefix_ + std::string( "surfaceBoundaryType" ) );
      CMBBoundaryType_                                     = parameters_.getParameter< uint_t >( prefix_ + std::string( "CMBBoundaryType" ) );
      temperatureExtrapolationOrder_                       = parameters_.getParameter< uint_t >( prefix_ + std::string( "temperatureExtrapolationOrder" ) );
      velocityExtrapolationOrder_                          = parameters_.getParameter< uint_t >( prefix_ + std::string( "velocityExtrapolationOrder" ) );
      BDFOrder_                                            = parameters_.getParameter< uint_t >( prefix_ + std::string( "BDFOrder" ) );
      viscosityDerivativeOrder_                            = parameters_.getParameter< uint_t >( prefix_ + std::string( "viscosityDerivativeOrder" ) );
      densityDerivativeOrder_                              = parameters_.getParameter< uint_t >( prefix_ + std::string( "densityDerivativeOrder" ) );
      maxSteps_                                            = parameters_.getParameter< uint_t >( prefix_ + std::string( "maxSteps" ) );
      maximumSaddlePointSolverIterations_                  = parameters_.getParameter< uint_t >( prefix_ + std::string( "maximumSaddlePointSolverIterations" ) );
      temperatureMassMaxIterations_                        = parameters_.getParameter< uint_t >( prefix_ + std::string( "temperatureMassMaxIterations" ) );
      loadCheckpointNumber_                                = parameters_.getParameter< uint_t >( prefix_ + std::string( "loadCheckpointNumber" ) );

      boundaryTolerance_                                   = parameters_.getParameter< real_t >( prefix_ + std::string( "boundaryTolerance" ) );
      plateVelocityScaling_                                = parameters_.getParameter< real_t >( prefix_ + std::string( "plateVelocityScaling" ) );
      minTimestepMyrs_                                     = parameters_.getParameter< real_t >( prefix_ + std::string( "minTimestepMyrs" ) );
      maxTimestepMyrs_                                     = parameters_.getParameter< real_t >( prefix_ + std::string( "maxTimestepMyrs" ) );
      CFL_                                                 = parameters_.getParameter< real_t >( prefix_ + std::string( "CFL" ) );
      relativeResidualToleranceOuterSaddlePointSolverLoop_ = parameters_.getParameter< real_t >( prefix_ + std::string( "relativeResidualToleranceOuterSaddlePointSolverLoop" ) );
      absoluteResidualToleranceOuterSaddlePointSolverLoop_ = parameters_.getParameter< real_t >( prefix_ + std::string( "absoluteResidualToleranceOuterSaddlePointSolverLoop" ) );
      temperatureMassRelativeTolerance_                    = parameters_.getParameter< real_t >( prefix_ + std::string( "temperatureMassRelativeTolerance" ) );
      temperatureMassAbsoluteTolerance_                    = parameters_.getParameter< real_t >( prefix_ + std::string( "temperatureMassAbsoluteTolerance" ) );
      SUPG_scaling_                                        = parameters_.getParameter< real_t >( prefix_ + std::string( "SUPG_scaling" ) );
      const_H_                                             = parameters_.getParameter< real_t >( prefix_ + std::string( "const_H" ) );
      const_alpha_                                         = parameters_.getParameter< real_t >( prefix_ + std::string( "const_alpha" ) );
      const_K_T_                                           = parameters_.getParameter< real_t >( prefix_ + std::string( "const_K_T" ) );
      const_k_                                             = parameters_.getParameter< real_t >( prefix_ + std::string( "const_k" ) );
      const_C_p_                                           = parameters_.getParameter< real_t >( prefix_ + std::string( "const_C_p" ) );

      usePlates_                                           = parameters_.getParameter< bool >( prefix_ + std::string( "usePlates" ) );
      loopPlateAge_                                        = parameters_.getParameter< bool >( prefix_ + std::string( "loopPlateAge" ) );
      randomInitialGuessU_                                 = parameters_.getParameter< bool >( prefix_ + std::string( "randomInitialGuessU" ) );
      randomInitialGuessP_                                 = parameters_.getParameter< bool >( prefix_ + std::string( "randomInitialGuessP" ) );      
      vtk_                                                 = parameters_.getParameter< bool >( prefix_ + std::string( "vtk" ) );     
      bp4_                                                 = parameters_.getParameter< bool >( prefix_ + std::string( "bp4" ) );     
      writeCheckpoints_                                    = parameters_.getParameter< bool >( prefix_ + std::string( "writeCheckpoints" ) ); 
      defaultToUsingAdiosCheckpoints_                      = parameters_.getParameter< bool >( prefix_ + std::string( "defaultToUsingAdiosCheckpoints" ) );  
      temperatureMassPrintInfo_                            = parameters_.getParameter< bool >( prefix_ + std::string( "temperatureMassPrintInfo" ) );    
      loadCheckpointOnStart_                               = parameters_.getParameter< bool >( prefix_ + std::string( "loadCheckpointOnStart" ) );   
      lowMemoryMode_                                       = parameters_.getParameter< bool >( prefix_ + std::string( "lowMemoryMode" ) );  
      alwaysDestroyTemporaryFunctions_                     = parameters_.getParameter< bool >( prefix_ + std::string( "alwaysDestroyTemporaryFunctions" ) );  
      useGlobalCFL_                                        = parameters_.getParameter< bool >( prefix_ + std::string( "useGlobalCFL" ) );  
      useMyrsInsteadOfTimeStepsAsOutputFrequency_          = parameters_.getParameter< bool >( prefix_ + std::string( "useMyrsInsteadOfTimeStepsAsOutputFrequency" ) ); 

      fileTopologies_                                      = parameters_.getParameter< std::string >( prefix_ + std::string( "fileTopologies" ) );
      fileReconstructions_                                 = parameters_.getParameter< std::string >( prefix_ + std::string( "fileReconstructions" ) );
      fileName_                                            = parameters_.getParameter< std::string >( prefix_ + std::string( "fileName" ) );
      // clang-format on

      // load parameter file to vector on root (used to store in checkpoints)
      WALBERLA_ROOT_SECTION()
      {
         if ( writeCheckpoints_ )
         {
            std::ifstream in( parameterfile.c_str(), std::ifstream::in );
            if ( !in.is_open() )
            {
               WALBERLA_ABORT( "Error opening parameter file." )
            }

            std::string line;

            while ( std::getline( in, line ) )
            {
               parameterFileVector_.push_back( line );
            }

            in.close();
         }
      }

      // set alwaysDestroyTempFunctions
      TempFunctionManager::instance()->setAlwaysDestroy( alwaysDestroyTemporaryFunctions_ );

      // create NondimensionalisationParameters
      ND_ = NondimensionalisationParameters( parameters_, prefix_ );

      // dependent parameters
      minDt_ = minTimestepMyrs_ * secondsPerMillionYears_ / ND_.tRef_;
      maxDt_ = maxTimestepMyrs_ * secondsPerMillionYears_ / ND_.tRef_;

      // define boundary types
      switch ( surfaceBoundaryType_ )
      {
      case 1: // 1 = Dirichlet
         SurfaceType_ = DirichletBoundary;
         break;
      case 2: // 2 = Neumann
         SurfaceType_ = NeumannBoundary;
         break;
      case 3: // 3 = FreeSlip
         SurfaceType_ = FreeslipBoundary;
         break;
      default: //default = Dirichlet
         SurfaceType_ = DirichletBoundary;
         break;
      }

      switch ( CMBBoundaryType_ )
      {
      case 1: // 1 = Dirichlet
         CMBType_ = DirichletBoundary;
         break;
      case 2: // 2 = Neumann
         CMBType_ = NeumannBoundary;
         break;
      case 3: // 3 = FreeSlip
         CMBType_ = FreeslipBoundary;
         break;
      default: //default = Dirichlet
         CMBType_ = DirichletBoundary;
         break;
      }

      // create storage, indicator functions and boundary conditions
      if ( dim_ == 2 )
      {
         nTan_ = parameters_.getParameter< uint_t >( prefix_ + std::string( "nTan2D" ) );
         nRad_ = parameters_.getParameter< uint_t >( prefix_ + std::string( "nRad2D" ) );

         if ( loadStorageFile_ )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "loading setupPrimitiveStorage file " << fileSetupPrimitiveStorage_ << "..." );
            storage_ = std::make_shared< PrimitiveStorage >( fileSetupPrimitiveStorage_, 1 );
            WALBERLA_LOG_INFO_ON_ROOT( "done..." );
         }
         else
         {
            storage_ = createAnnulusStorage( nTan_, nRad_, ND_.radiusCMB_, ND_.radiusSurface_, blending_, boundaryTolerance_ );
         }

         surfaceFct_ = createAnnulusSurfaceBoundaryFct( ND_.radiusSurface_, boundaryTolerance_ );
         CMBFct_     = createAnnulusCMBBoundaryFct( ND_.radiusCMB_, boundaryTolerance_ );

         BC_ = createAnnulusBoundaryConditions( SurfaceType_, CMBType_ );
      }
      else if ( dim_ == 3 )
      {
         nTan_ = parameters_.getParameter< uint_t >( prefix_ + std::string( "nTan3D" ) );
         nRad_ = parameters_.getParameter< uint_t >( prefix_ + std::string( "nRad3D" ) );

         if ( loadStorageFile_ )
         {
            WALBERLA_LOG_INFO_ON_ROOT( "loading setupPrimitiveStorage file " << fileSetupPrimitiveStorage_ << "..." );
            storage_ = std::make_shared< PrimitiveStorage >( fileSetupPrimitiveStorage_, 1 );
            WALBERLA_LOG_INFO_ON_ROOT( "done..." );
         }
         else
         {
            storage_ =
                createSphericalShellStorage( nTan_, nRad_, ND_.radiusCMB_, ND_.radiusSurface_, blending_, boundaryTolerance_ );
         }

         surfaceFct_ = createSphericalShellSurfaceBoundaryFct( ND_.radiusSurface_, boundaryTolerance_ );
         CMBFct_     = createSphericalShellCMBBoundaryFct( ND_.radiusCMB_, boundaryTolerance_ );

         BC_ = createSphericalShellBoundaryConditions( SurfaceType_, CMBType_ );
      }
      else
      {
         WALBERLA_ABORT( "Dimension " << dim_ << " not supported!" );
      }

      // check BDF, MMOC combination
      if ( BDFOrder_ == 0 )
      {
         WALBERLA_ABORT( "BDFOrder " << BDFOrder_ << " not supported!" );
      }
      if ( BDFOrder_ > 2 && MMOC_ )
      {
         WALBERLA_ABORT( "BDFOrder " << BDFOrder_ << " not supported with MMOC!" );
      }

      // storage properties
      nDoFsAdvectionDiffusionGlobal_ =
          hyteg::numberOfGlobalDoFs< hyteg::FunctionTrait< P2Function< real_t > >::Tag >( *storage_, maxLevel_ );
      nDoFsSaddlePointGlobal_ =
          hyteg::numberOfGlobalDoFs< hyteg::FunctionTrait< P2P1TaylorHoodFunction< real_t > >::Tag >( *storage_, maxLevel_ );
      nDoFsAdvectionDiffusionLocal_ =
          hyteg::numberOfLocalDoFs< hyteg::FunctionTrait< P2Function< real_t > >::Tag >( *storage_, maxLevel_ );
      nDoFsSaddlePointLocal_ =
          hyteg::numberOfLocalDoFs< hyteg::FunctionTrait< P2P1TaylorHoodFunction< real_t > >::Tag >( *storage_, maxLevel_ );
      nDoFsAdvectionDiffusionLocalMax_ = walberla::mpi::allReduce( nDoFsAdvectionDiffusionLocal_, walberla::mpi::MAX );
      nDoFsSaddlePointLocalMax_        = walberla::mpi::allReduce( nDoFsSaddlePointLocal_, walberla::mpi::MAX );

      hMin_ = hyteg::MeshQuality::getMinimalEdgeLength( storage_, maxLevel_ );
      hMax_ = hyteg::MeshQuality::getMaximalEdgeLength( storage_, maxLevel_ );
      if ( dim == 3 )
      {
         nMacroElementsGlobal_ = storage_->getNumberOfGlobalCells();
         nMacroElementsLocal_  = storage_->getNumberOfLocalCells();
      }
      else
      {
         nMacroElementsGlobal_ = storage_->getNumberOfGlobalFaces();
         nMacroElementsLocal_  = storage_->getNumberOfLocalCells();
      }
      nMacroElementsLocalMax_ = walberla::mpi::allReduce( nMacroElementsLocal_, walberla::mpi::MAX );

      nDoFsAdvectionDiffusionLocalAvg_ =
          real_c( walberla::mpi::allReduce( nDoFsAdvectionDiffusionLocal_, walberla::mpi::SUM ) ) / real_c( size_ );
      nDoFsSaddlePointLocalAvg_ =
          real_c( walberla::mpi::allReduce( nDoFsSaddlePointLocal_, walberla::mpi::SUM ) ) / real_c( size_ );
      real_t localVarAdvectionDiffusion =
          std::pow( real_c( nDoFsAdvectionDiffusionLocal_ ) - nDoFsAdvectionDiffusionLocalAvg_, 2 );
      real_t localVarSaddlePoint = std::pow( real_c( nDoFsSaddlePointLocal_ ) - nDoFsSaddlePointLocalAvg_, 2 );

      nDoFsAdvectionDiffusionLocalVar_ =
          walberla::mpi::allReduce( localVarAdvectionDiffusion, walberla::mpi::SUM ) / real_c( size_ );
      nDoFsSaddlePointLocalVar_ = walberla::mpi::allReduce( localVarSaddlePoint, walberla::mpi::SUM ) / real_c( size_ );

      // define freeslip projection
      if ( dim == 3 )
      {
         surfaceNormal_ = [&]( const Point3D& x, Point3D& normal ) {
            real_t radius = x.norm();
            if ( std::abs( radius - ND_.radiusSurface_ ) < std::abs( radius - ND_.radiusCMB_ ) )
            {
               normal = Point3D( { x[0] / radius, x[1] / radius, x[2] / radius } );
            }
            else
            {
               normal = Point3D( { -x[0] / radius, -x[1] / radius, -x[2] / radius } );
            }
         };
      }
      else
      {
         surfaceNormal_ = [&]( const Point3D& x, Point3D& normal ) {
            real_t radius = x.norm();
            if ( std::abs( radius - ND_.radiusSurface_ ) < std::abs( radius - ND_.radiusCMB_ ) )
            {
               normal = Point3D( { x[0] / radius, x[1] / radius, real_c( 0 ) } );
            }
            else
            {
               normal = Point3D( { -x[0] / radius, -x[1] / radius, real_c( 0 ) } );
            }
         };
      }

      projection_ = std::make_shared< P2ProjectNormalOperator >( storage_, minLevel_, maxLevel_, surfaceNormal_ );

      // define other std::functions
      zeroFct_ = []( const Point3D& x ) {
         WALBERLA_UNUSED( x );
         return real_c( 0 );
      };
      randFunc_ = []( const Point3D& ) { return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) ); };

      // setup FEM functions
      up_extra_ = createSharedP2P1TaylorHoodFunction( "up_extra", storage_, minLevel_, maxLevel_, BC_ );
      fg_       = createSharedP2P1TaylorHoodFunction( "fg", storage_, minLevel_, maxLevel_, BC_ );
      if ( !lowMemoryMode_ )
      {
         tmpSaddlePoint_ = createSharedP2P1TaylorHoodFunction( "tmpSaddlePoint", storage_, minLevel_, maxLevel_, BC_ );
      }

      T_d_     = createSharedP2TemperatureFunction( "T_d", storage_, minLevel_, maxLevel_, BC_ );
      T_extra_ = createSharedP2TemperatureFunction( "T_extra", storage_, minLevel_, maxLevel_, BC_ );
      h_       = createSharedP2TemperatureFunction( "h", storage_, minLevel_, maxLevel_, BC_ );
      if ( !lowMemoryMode_ )
      {
         tmpTemp_ = createSharedP2TemperatureFunction( "tmpTemp", storage_, minLevel_, maxLevel_, BC_ );
         if ( BDFOrder_ >= 2 )
         {
            tmpTemp2_          = createSharedP2TemperatureFunction( "tmpTemp2", storage_, minLevel_, maxLevel_, BC_ );
            mmocDummyFunction_ = createSharedP2TemperatureFunction( "mmocDummyFunction_", storage_, minLevel_, maxLevel_, BC_ );
         }
      }

      eta_extra_ = std::make_shared< ViscosityFunctionType_ >(
          "eta extra", storage_, minLevel_, maxLevel_, BoundaryCondition::createAllInnerBC() );
      rho_extra_ = std::make_shared< DensityFunctionType_ >(
          "rho extra", storage_, minLevel_, maxLevel_, BoundaryCondition::createAllInnerBC() );

      if constexpr ( CalculateInvRhoAndEta_ )
      {
         inv_eta_ = std::make_shared< ViscosityFunctionType_ >(
             "inv eta", storage_, minLevel_, maxLevel_, BoundaryCondition::createAllInnerBC() );
         inv_rho_ = std::make_shared< DensityFunctionType_ >(
             "inv rho", storage_, minLevel_, maxLevel_, BoundaryCondition::createAllInnerBC() );
         inv_eta_extra_ = std::make_shared< ViscosityFunctionType_ >(
             "inv eta extra", storage_, minLevel_, maxLevel_, BoundaryCondition::createAllInnerBC() );
         inv_rho_extra_ = std::make_shared< DensityFunctionType_ >(
             "inv rho extra", storage_, minLevel_, maxLevel_, BoundaryCondition::createAllInnerBC() );
      }

      // setup function histories
      up_ = std::make_shared< FunctionHistory< P2P1TaylorHoodFunction< real_t >, real_t > >(
          std::max( velocityExtrapolationOrder_ + 1, BDFOrder_ ), minLevel_, maxLevel_ );
      for ( uint_t i = 0; i < up_->getMemoryCapacity(); i++ )
      {
         std::stringstream fName;
         fName << "up history " << i;
         auto upFct = createSharedP2P1TaylorHoodFunction( fName.str(), storage_, minLevel_, maxLevel_, BC_ );
         up_->addFunction( upFct );
      }

      T_ = std::make_shared< FunctionHistory< P2Function< real_t >, real_t > >(
          std::max( temperatureExtrapolationOrder_, BDFOrder_ ) + 1, minLevel_, maxLevel_ );
      for ( uint_t i = 0; i < T_->getMemoryCapacity(); i++ )
      {
         std::stringstream fName;
         fName << "T history " << i;
         auto TFct = createSharedP2TemperatureFunction( fName.str(), storage_, minLevel_, maxLevel_, BC_ );
         T_->addFunction( TFct );
      }

      rho_ = std::make_shared< FunctionHistory< DensityFunctionType_, real_t > >(
          densityDerivativeOrder_ + 1, minLevel_, maxLevel_ );
      for ( uint_t i = 0; i < rho_->getMemoryCapacity(); i++ )
      {
         std::stringstream fName;
         fName << "rho history " << i;
         auto rhoFct = std::make_shared< DensityFunctionType_ >(
             fName.str(), storage_, minLevel_, maxLevel_, BoundaryCondition::createAllInnerBC() );
         rho_->addFunction( rhoFct );
      }

      eta_ = std::make_shared< FunctionHistory< ViscosityFunctionType_, real_t > >(
          viscosityDerivativeOrder_ + 1, minLevel_, maxLevel_ );
      for ( uint_t i = 0; i < eta_->getMemoryCapacity(); i++ )
      {
         std::stringstream fName;
         fName << "eta history " << i;
         auto etaFct = std::make_shared< ViscosityFunctionType_ >(
             fName.str(), storage_, minLevel_, maxLevel_, BoundaryCondition::createAllInnerBC() );
         eta_->addFunction( etaFct );
      }

      // init mmoc
      if ( MMOC_ )
      {
         if ( !lowMemoryMode_ )
         {
            mmocTransport_ = std::make_shared< MMOCTransport< P2Function< real_t > > >(
                storage_, minLevel_, maxLevel_, TimeSteppingScheme::RK4, lowMemoryMode_ );
         }
      }

      // init plate oracle
      if ( usePlates_ )
      {
         plateOracle_  = std::make_shared< terraneo::plates::PlateVelocityProvider >( fileTopologies_, fileReconstructions_ );
         real_t minAge = plateOracle_->getMinAge();
         real_t maxAge = plateOracle_->getMaxAge();
         WALBERLA_LOG_INFO_ON_ROOT( "Minage: " << minAge << " ; Maxage: " << maxAge );

         if ( loopPlateAge_ && getCurrentPlateAge() >= plateOracle_->getMaxAge() )
         {
            currentPlateTime_ = currentPlateTime_ - plateOracle_->getMaxAge() * secondsPerMillionYears_ / ND_.tRef_;
         }
      }

      // init temperature mass
      temperatureMass_ = std::make_shared< TemperatureMassOperatorType_ >( storage_, minLevel_, maxLevel_ );
      std::function< real_t( const Point3D& ) > P2P1MassScaling = [=]( const hyteg::Point3D& x ) {
         WALBERLA_UNUSED( x );
         return real_c( 1.0 );
      };
      hyteg::forms::p2_to_p1_k_mass_blending_q5 FormP2P1Mass( P2P1MassScaling );
      temperatureP2ToP1Mass_ =
          std::make_shared< P2ToP1ElementwiseBlendingKMassOperator >( storage_, minLevel_, maxLevel_, FormP2P1Mass );
      temperatureRestriction_ = std::make_shared< TemperatureRestrictionOperatorType_ >();
      temperatureMassSolver_  = std::make_shared< hyteg::CGSolver< TemperatureMassOperatorType_ > >(
          storage_,
          minLevel_,
          maxLevel_,
          temperatureMassMaxIterations_,
          temperatureMassRelativeTolerance_,
          temperatureMassAbsoluteTolerance_,
          std::make_shared< hyteg::IdentityPreconditioner< TemperatureMassOperatorType_ > >(),
          lowMemoryMode_ );
      temperatureMassSolver_->setName( "CG Temperature Mass" );
      temperatureMassSolver_->setPrintInfo( temperatureMassPrintInfo_ );

      // init operator updater
      operatorUpdater_ = std::make_shared< IdentityOperatorUpdater >();

      // checkpoint
      checkpointNotLoaded_ = true;

      // init vtk output
      if ( vtk_ )
      {
         if ( useMyrsInsteadOfTimeStepsAsOutputFrequency_ )
         {
            vtkOutput_ = std::make_shared< VTKOutput >( "./vtk", fileName_, storage_, 1 );
         }
         else
         {
            vtkOutput_ = std::make_shared< VTKOutput >( "./vtk", fileName_, storage_, writeFrequencyOutput_ );
         }

         for ( uint_t dimension = 0; dimension < ( storage_->hasGlobalCells() ? 3 : 2 ); dimension++ )
         {
            vtkOutput_->add( up_->getFunctionByIndex( 0 ).uvw()[dimension] );
            vtkOutput_->add( up_extra_->uvw()[dimension] );
            vtkOutput_->add( fg_->uvw()[dimension] );
         }
         vtkOutput_->add( up_->getFunctionByIndex( 0 ).p() );
         vtkOutput_->add( up_extra_->p() );
         vtkOutput_->add( fg_->p() );
         vtkOutput_->add( *T_d_ );
         vtkOutput_->add( T_->getFunctionByIndex( 0 ) );
         vtkOutput_->add( *T_extra_ );
         vtkOutput_->add( *h_ );
         vtkOutput_->add( rho_->getFunctionByIndex( 0 ) );
         vtkOutput_->add( eta_->getFunctionByIndex( 0 ) );
         vtkOutput_->add( *eta_extra_ );
         vtkOutput_->add( *rho_extra_ );
      }

#ifdef HYTEG_BUILD_WITH_ADIOS2
      {
         if ( bp4_ )
         {
            bp4Output_ = std::make_shared< AdiosWriter >( "./vtk", fileName_, storage_ );

            for ( uint_t dimension = 0; dimension < ( storage_->hasGlobalCells() ? 3 : 2 ); dimension++ )
            {
               bp4Output_->add( up_->getFunctionByIndex( 0 ).uvw()[dimension] );
               bp4Output_->add( up_extra_->uvw()[dimension] );
               bp4Output_->add( fg_->uvw()[dimension] );
            }
            bp4Output_->add( up_->getFunctionByIndex( 0 ).p() );
            bp4Output_->add( up_extra_->p() );
            bp4Output_->add( fg_->p() );
            bp4Output_->add( *T_d_ );
            bp4Output_->add( T_->getFunctionByIndex( 0 ) );
            bp4Output_->add( *T_extra_ );
            bp4Output_->add( *h_ );
            bp4Output_->add( rho_->getFunctionByIndex( 0 ) );
            bp4Output_->add( eta_->getFunctionByIndex( 0 ) );
            bp4Output_->add( *eta_extra_ );
            bp4Output_->add( *rho_extra_ );
         }
      }
#endif
   };

   // to be called after setting the temperature, density and viscosity model
   // but before constructing and setting your operators and solvers
   void init()
   {
      // initial guess for uvw
      up_->newState( real_c( dt_ ) );
      if ( randomInitialGuessU_ )
      {
         up_->getState( 0 ).uvw().interpolate( randFunc_, maxLevel_, solverFlag_ );
      }
      else
      {
         up_->getState( 0 ).uvw().interpolate( real_c( 0 ), maxLevel_, solverFlag_ );
      }
      up_->getState( 0 ).uvw().interpolate( real_c( 0 ), maxLevel_, dirichletFlag_ );

      // initial guess for p
      if ( randomInitialGuessP_ )
      {
         up_->getState( 0 ).p().interpolate( randFunc_, maxLevel_, solverFlag_ );
      }
      else
      {
         up_->getState( 0 ).p().interpolate( real_c( 0 ), maxLevel_, solverFlag_ );
      }
      up_->getState( 0 ).p().interpolate( real_c( 0 ), maxLevel_, dirichletFlag_ );

      // make sure the new state / initial guess obeys the boundary conditions before solving
      enforceSaddlePointBoundaryConditions( up_->getState( 0 ), maxLevel_ );

      // initial state for T_;
      T_->newState( real_c( dt_ ) );
      initialTemperatureModel_->interpolate( T_->getState( 0 ), maxLevel_, allFlag_ );

      // propagate temperature to coarser levels
      propagateTemperatureToCoarserLevel( T_->getState( 0 ), maxLevel_, minLevel_ );

      // initial extrapolations, do an incompressible first stokes solve
      up_extra_->interpolate( real_c( 0 ), maxLevel_, solverFlag_ );
      T_extra_->assign( { real_c( 1 ) }, { T_->getState( 0 ) }, maxLevel_, allFlag_ );

      // init dynamic temperature
      referenceTemperature_ = referenceTemperatureModel_->getTemperatureFctNoCheck();
      dynamicTemperature_   = [&]( const hyteg::Point3D& x, const std::vector< real_t >& fields ) {
         return fields[0] - referenceTemperature_( x );
      };

      // initial interpolate density
      rho_->newState( 0 );
      if constexpr ( InterpolateDensity_ )
      {
         for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
         {
            interpolateDensity( rho_->getState( 0 ), T_->getState( 0 ), level, allFlag_ );
            rho_extra_->assign( { real_c( 1 ) }, { rho_->getState( 0 ) }, level, allFlag_ );

            if constexpr ( CalculateInvRhoAndEta_ )
            {
               inv_rho_->assign( { real_c( 1 ) }, { rho_->getState( 0 ) }, level, allFlag_ );
               inv_rho_->invertElementwise( level, allFlag_, false );

               inv_rho_extra_->assign( { real_c( 1 ) }, { *rho_extra_ }, level, allFlag_ );
               inv_rho_extra_->invertElementwise( level, allFlag_, false );
            }
         }
      }

      // initial interpolate viscosity
      eta_->newState( 0 );
      if constexpr ( InterpolateViscosity_ )
      {
         for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
         {
            interpolateViscosity( eta_->getState( 0 ), T_->getState( 0 ), level, allFlag_ );
            eta_extra_->assign( { real_c( 1 ) }, { eta_->getState( 0 ) }, level, allFlag_ );

            if constexpr ( CalculateInvRhoAndEta_ )
            {
               inv_eta_->assign( { real_c( 1 ) }, { eta_->getState( 0 ) }, level, allFlag_ );
               inv_eta_->invertElementwise( level, allFlag_, false );

               inv_eta_extra_->assign( { real_c( 1 ) }, { *eta_extra_ }, level, allFlag_ );
               inv_eta_extra_->invertElementwise( level, allFlag_, false );
            }
         }
      }
   }

   void step()
   {
      if ( loadCheckpointOnStart_ && checkpointNotLoaded_ )
      {
         loadCheckpoint( fileName_, loadCheckpointNumber_ );
         checkpointNotLoaded_ = false;

         // update operators after checkpoint load
         operatorUpdater_->updateAfterCheckPointLoad();
      }
      else
      {
         // ############################
         // #### Saddle point solve ####
         // ############################

         // use extrapolation as initial guess
         // up_extra_ already fulfills the boundary conditions
         if ( stepCounter_ > 0 )
         {
            up_->getState( 0 ).assign( { real_c( 1 ) }, { *up_extra_ }, maxLevel_, solverFlag_ );
         }

         // calculate saddle point RHS
         calculateSaddlePointRHS();

         // update operators before saddle point solve
         operatorUpdater_->updateBeforeSaddlePointSolve();

         // saddle point solver loop
         hyteg::Point3D initResidual = calculateSaddlePointResidual( up_->getState( 0 ), *fg_, maxLevel_ );
         WALBERLA_LOG_INFO_ON_ROOT( "Initial residual: " << initResidual[0] << " ; Initial residual U: " << initResidual[1]
                                                         << " ; Initial residual P: " << initResidual[2] );

         hyteg::Point3D prevResidual = initResidual;

         walberla::WcTimer saddlePointTimer_;

         for ( uint_t i = 0; i <= maximumSaddlePointSolverIterations_; i++ )
         {
            // check residual
            hyteg::Point3D residual    = calculateSaddlePointResidual( up_->getState( 0 ), *fg_, maxLevel_ );
            hyteg::Point3D relResidual = {
                residual[0] / initResidual[0], residual[1] / initResidual[1], residual[2] / initResidual[2] };
            hyteg::Point3D reductionRate = {
                residual[0] / prevResidual[0], residual[1] / prevResidual[1], residual[2] / prevResidual[2] };
            prevResidual = residual;

            // output
            WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
                "[Iteration Callback] iter %3d | abs. res: %10.5e | abs. res U: %10.5e | abs. res P: %10.5e | rel. res: %10.5e | rel. res U: %10.5e | rel. res P: %10.5e | conv. rate: %10.5e",
                i,
                residual[0],
                residual[1],
                residual[2],
                relResidual[0],
                relResidual[1],
                relResidual[2],
                reductionRate[0] ) );

            if ( relResidual[0] < relativeResidualToleranceOuterSaddlePointSolverLoop_ )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "Relative residual tolerance of "
                                          << std::scientific << relativeResidualToleranceOuterSaddlePointSolverLoop_
                                          << " reached after " << i << " iterations. " );
               break;
            }

            if ( residual[0] < absoluteResidualToleranceOuterSaddlePointSolverLoop_ )
            {
               WALBERLA_LOG_INFO_ON_ROOT( "Absolute residual tolerance of "
                                          << std::scientific << absoluteResidualToleranceOuterSaddlePointSolverLoop_
                                          << " reached after " << i << " iterations. " );
               break;
            }

            saddlePointTimer_.start();
            saddlePointSolver_->solve( *saddlePointOperator_, up_->getState( 0 ), *fg_, maxLevel_ );
            saddlePointTimer_.end();
         }

         WALBERLA_LOG_INFO_ON_ROOT( "Total Saddlepoint solve time: " << saddlePointTimer_.total() << " seconds." );

         // update operators after saddle point solve
         operatorUpdater_->updateAfterSaddlePointSolve();

         // ################
         // #### Output ####
         // ################
         const bool currentTimeZero = std::fpclassify( currentTime_ ) == FP_ZERO;

         if ( vtk_ )
         {
            if ( useMyrsInsteadOfTimeStepsAsOutputFrequency_ )
            {
               if ( saveNextOutput_ || currentTimeZero )
               {
                  vtkOutput_->write( maxLevel_, stepCounter_ );
               }
            }
            else
            {
               vtkOutput_->write( maxLevel_, stepCounter_ );
            }
         }
#ifdef HYTEG_BUILD_WITH_ADIOS2
         {
            if ( bp4_ )
            {
               if ( useMyrsInsteadOfTimeStepsAsOutputFrequency_ )
               {
                  if ( saveNextOutput_ || currentTimeZero )
                  {
                     bp4Output_->write( maxLevel_, stepCounter_ );
                  }
               }
               else
               {
                  if ( writeFrequencyOutput_ > 0 && stepCounter_ % writeFrequencyOutput_ == 0 )
                  {
                     bp4Output_->write( maxLevel_, stepCounter_ );
                  }
               }
            }
         }
#endif

         // ##########################
         // #### Write Checkpoint ####
         // ##########################
         if ( writeCheckpoints_ )
         {
            if ( useMyrsInsteadOfTimeStepsAsOutputFrequency_ )
            {
               if ( saveNextCheckpoint_ || currentTimeZero )
               {
                  createCheckpoint( fileName_ );
               }
            }
            else
            {
               if ( writeFrequencyCheckpoint_ > 0 && stepCounter_ % writeFrequencyCheckpoint_ == 0 )
               {
                  createCheckpoint( fileName_ );
               }
            }
         }
      }

      // ###################################
      // #### Advection diffusion solve ####
      // ###################################

      WALBERLA_LOG_INFO_ON_ROOT( "---------------- Step " << stepCounter_ + 1 << " ----------------" );

      // calculate new time step via CFL
      if ( useGlobalCFL_ )
      {
         std::shared_ptr< P2P1TaylorHoodFunction< real_t > > tmpSaddlePointStep;

         if ( !lowMemoryMode_ )
         {
            tmpSaddlePointStep = tmpSaddlePoint_;
         }
         else
         {
            tmpSaddlePointStep = getTemporaryFunction< P2P1TaylorHoodFunction< real_t > >( storage_, minLevel_, maxLevel_ );

            tmpSaddlePointStep->uvw().component( 0 ).setBoundaryCondition( BC_.bcVelocityX_ );
            tmpSaddlePointStep->uvw().component( 1 ).setBoundaryCondition( BC_.bcVelocityY_ );
            if ( storage_->hasGlobalCells() )
            {
               tmpSaddlePointStep->uvw().component( 2 ).setBoundaryCondition( BC_.bcVelocityZ_ );
            }
         }

         dt_ = std::max( std::min( CFLTimestep( up_->getState( 0 ).uvw(), maxLevel_, tmpSaddlePointStep->uvw(), true ), maxDt_ ),
                         minDt_ );

         WALBERLA_LOG_INFO_ON_ROOT( "New time step (global estimation): " << dt_ );
      }
      else
      {
         dt_ = std::max( std::min( CFLTimestep( up_->getState( 0 ).uvw(), maxLevel_, up_->getState( 0 ).uvw(), false ), maxDt_ ),
                         minDt_ );

         WALBERLA_LOG_INFO_ON_ROOT( "New time step (local estimation): " << dt_ );
      }

      // update age
      real_t stepSize = dt_ * ND_.tRef_ / secondsPerMillionYears_;
      real_t lastAge  = getCurrentAge();

      currentTime_ += dt_;

      real_t currAge = getCurrentAge();

      if ( useMyrsInsteadOfTimeStepsAsOutputFrequency_ )
      {
         if ( std::floor( currAge / real_c( writeFrequencyOutput_ ) ) > std::floor( lastAge / real_c( writeFrequencyOutput_ ) ) )
         {
            saveNextOutput_ = true;
         }
         else
         {
            saveNextOutput_ = false;
         }

         if ( std::floor( currAge / real_c( writeFrequencyCheckpoint_ ) ) >
              std::floor( lastAge / real_c( writeFrequencyCheckpoint_ ) ) )
         {
            saveNextCheckpoint_ = true;
         }
         else
         {
            saveNextCheckpoint_ = false;
         }
      }

      if ( usePlates_ )
      {
         currentPlateTime_ += dt_ / plateVelocityScaling_;

         if ( loopPlateAge_ && getCurrentPlateAge() >= plateOracle_->getMaxAge() )
         {
            currentPlateTime_ = currentPlateTime_ - plateOracle_->getMaxAge() * secondsPerMillionYears_ / ND_.tRef_;
         }
      }

      real_t currPlateAge = getCurrentPlateAge();

      WALBERLA_LOG_INFO_ON_ROOT( "dt_: " << dt_ << " ; stepSize: " << stepSize << " Myrs ; currentAge: " << currAge
                                         << " Myrs ; currentPlateAge: " << currPlateAge << " Myrs" );

      // calculate extrapolations
      if constexpr ( timeDependent_ )
      {
         up_->extrapolate( velocityExtrapolationOrder_, dt_, *up_extra_, maxLevel_, solverFlag_ );
         T_->extrapolate( temperatureExtrapolationOrder_, dt_, *T_extra_, maxLevel_, solverFlag_ );
      }
      else
      {
         up_extra_->assign( { real_c( 1 ) }, { up_->getState( 0 ) }, maxLevel_, allFlag_ );
         T_extra_->assign( { real_c( 1 ) }, { T_->getState( 0 ) }, maxLevel_, allFlag_ );
      }

      // make sure the extrapolations obey the boundary conditions
      // since we already updated the currentPlateTime_ this also sets the
      // boundary conditions of the velocity extrapolation to the correct
      // plate velocities in the future
      enforceSaddlePointBoundaryConditions( *up_extra_, maxLevel_ );
      enforceTemperatureBoundaryConditions( *T_extra_, maxLevel_ );

      // update operators after up_ extrapolation recalculation
      operatorUpdater_->updateAfterUpExtrapolationRecalculation();

      // propagate temperature extrapolation to coarser levels
      propagateTemperatureToCoarserLevel( *T_extra_, maxLevel_, minLevel_ );

      // update operators after temperature extrapolation recalculation
      operatorUpdater_->updateAfterTemperatureExtrapolationRecalculation();

      // update density extrapolation
      if constexpr ( InterpolateDensity_ )
      {
         for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
         {
            interpolateDensity( *rho_extra_, *T_extra_, level, allFlag_ );

            if constexpr ( CalculateInvRhoAndEta_ )
            {
               inv_rho_extra_->assign( { real_c( 1 ) }, { *rho_extra_ }, level, allFlag_ );
               inv_rho_extra_->invertElementwise( level, allFlag_, false );
            }
         }
      }

      // update operators after density extrapolation recalculation
      operatorUpdater_->updateAfterDensityExtrapolationRecalculation();

      // update viscosity extrapolation
      if constexpr ( InterpolateViscosity_ )
      {
         for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
         {
            interpolateViscosity( *eta_extra_, *T_extra_, level, allFlag_ );

            if constexpr ( CalculateInvRhoAndEta_ )
            {
               inv_eta_extra_->assign( { real_c( 1 ) }, { *eta_extra_ }, level, allFlag_ );
               inv_eta_extra_->invertElementwise( level, allFlag_, false );
            }
         }
      }

      // update operators after viscosity extrapolation recalculation
      operatorUpdater_->updateAfterViscosityExtrapolationRecalculation();

      // new T state
      T_->newState( dt_ );

      // mmoc transport
      std::shared_ptr< P2Function< real_t > > tmpTempStep;
      std::shared_ptr< P2Function< real_t > > tmpTempStep2;
      std::shared_ptr< P2Function< real_t > > mmocDummy;
      if constexpr ( MMOC_ )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "Perform MMOC advection step..." );

         if ( !lowMemoryMode_ )
         {
            tmpTempStep = tmpTemp_;
            if ( BDFOrder_ >= 2 )
            {
               tmpTempStep2 = tmpTemp2_;
               mmocDummy    = mmocDummyFunction_;
            }
         }
         else
         {
            tmpTempStep = getTemporaryFunction< P2Function< real_t > >( storage_, minLevel_, maxLevel_ );
            tmpTempStep->setBoundaryCondition( BC_.bcTemperature_ );

            if ( BDFOrder_ >= 2 && T_->stateOffsetAvailable( 2 ) )
            {
               tmpTempStep2 = getTemporaryFunction< P2Function< real_t > >( storage_, minLevel_, maxLevel_ );
               tmpTempStep2->setBoundaryCondition( BC_.bcTemperature_ );

               mmocDummy = getTemporaryFunction< P2Function< real_t > >( storage_, minLevel_, maxLevel_ );
               mmocDummy->setBoundaryCondition( BC_.bcTemperature_ );
            }
         }

         // save state 1
         tmpTempStep->assign( { real_c( 1 ) }, { T_->getState( 1 ) }, maxLevel_, allFlag_ );
         // save state 2
         if ( BDFOrder_ >= 2 && T_->stateOffsetAvailable( 2 ) )
         {
            tmpTempStep2->assign( { real_c( 1 ) }, { T_->getState( 2 ) }, maxLevel_, allFlag_ );
         }

         std::shared_ptr< MMOCTransport< P2Function< real_t > > > localMMOC;
         if ( !lowMemoryMode_ )
         {
            localMMOC = mmocTransport_;
         }
         else
         {
            localMMOC = std::make_shared< MMOCTransport< P2Function< real_t > > >(
                storage_, minLevel_, maxLevel_, TimeSteppingScheme::RK4, lowMemoryMode_ );
         }

         localMMOC->step( T_->getState( 1 ),
                          up_extra_->uvw(),
                          up_->getState( 0 ).uvw(),
                          maxLevel_,
                          allFlag_,
                          dt_,
                          1,
                          true,
                          true,
                          false,
                          true );
         enforceTemperatureBoundaryConditions( T_->getState( 1 ), maxLevel_ );

         if ( BDFOrder_ >= 2 && T_->stateOffsetAvailable( 2 ) )
         {
            localMMOC->step(
                *mmocDummy, up_extra_->uvw(), up_->getState( 0 ).uvw(), maxLevel_, allFlag_, dt_, 1, true, true, false, true );
            localMMOC->step( T_->getState( 2 ),
                             up_->getState( 0 ).uvw(),
                             up_->getState( 1 ).uvw(),
                             maxLevel_,
                             allFlag_,
                             up_->getStateStepSize( 0 ),
                             1,
                             false,
                             true,
                             false,
                             true );
            enforceTemperatureBoundaryConditions( T_->getState( 2 ), maxLevel_ );
         }
      }

      // calculate RHS
      calculateAdvectionDiffusionRHS();

      if constexpr ( MMOC_ )
      {
         // restore old state 1 in case of mmoc
         T_->getState( 1 ).assign( { real_c( 1 ) }, { *tmpTempStep }, maxLevel_, allFlag_ );
         // restore old state 2 in case of mmoc
         if ( BDFOrder_ >= 2 && T_->stateOffsetAvailable( 2 ) )
         {
            T_->getState( 2 ).assign( { real_c( 1 ) }, { *tmpTempStep2 }, maxLevel_, allFlag_ );
         }
      }

      // use extrapolation as initial guess
      // T_extra_ already fulfills the boundary conditions
      T_->getState( 0 ).assign( { real_c( 1 ) }, { *T_extra_ }, maxLevel_, allFlag_ );

      // solve advection diffusion
      walberla::WcTimer advectionDiffusionTimer_;
      advectionDiffusionTimer_.start();
      advectionDiffusionSolver_->solve( *advectionDiffusionOperator_, T_->getState( 0 ), *h_, maxLevel_ );
      advectionDiffusionTimer_.end();

      WALBERLA_LOG_INFO_ON_ROOT( "Total Advection Diffusion solve time: " << advectionDiffusionTimer_.total() << " seconds." );

      // propagate temperature to coarser levels
      propagateTemperatureToCoarserLevel( T_->getState( 0 ), maxLevel_, minLevel_ );

      // update operators after advection diffusion solve
      operatorUpdater_->updateAfterAdvectionDiffusionSolve();

      // new up state
      up_->newState( dt_ );

      // interpolate density
      rho_->newState( dt_ );
      if constexpr ( InterpolateDensity_ )
      {
         for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
         {
            interpolateDensity( rho_->getState( 0 ), T_->getState( 0 ), level, allFlag_ );

            if constexpr ( CalculateInvRhoAndEta_ )
            {
               inv_rho_->assign( { real_c( 1 ) }, { rho_->getState( 0 ) }, level, allFlag_ );
               inv_rho_->invertElementwise( level, allFlag_, false );
            }
         }
      }

      // update operators after density recalculation
      operatorUpdater_->updateAfterDensityRecalculation();

      // interpolate viscosity
      eta_->newState( dt_ );
      if constexpr ( InterpolateViscosity_ )
      {
         for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
         {
            interpolateViscosity( eta_->getState( 0 ), T_->getState( 0 ), level, allFlag_ );

            if constexpr ( CalculateInvRhoAndEta_ )
            {
               inv_eta_->assign( { real_c( 1 ) }, { eta_->getState( 0 ) }, level, allFlag_ );
               inv_eta_->invertElementwise( level, allFlag_, false );
            }
         }
      }

      // update operators after viscosity recalculation
      operatorUpdater_->updateAfterViscosityRecalculation();

      // new step
      stepCounter_++;
   }

   void runSimulation()
   {
      WALBERLA_LOG_INFO_ON_ROOT( "---------------- Step " << stepCounter_ << " ----------------" );
      while ( stepCounter_ < maxSteps_ )
      {
         step();
      }

      // ################
      // #### Output ####
      // ################
      const bool currentTimeZero = std::fpclassify( currentTime_ ) == FP_ZERO;

      if ( vtk_ )
      {
         if ( useMyrsInsteadOfTimeStepsAsOutputFrequency_ )
         {
            if ( saveNextOutput_ || currentTimeZero )
            {
               vtkOutput_->write( maxLevel_, stepCounter_ );
            }
         }
         else
         {
            vtkOutput_->write( maxLevel_, stepCounter_ );
         }
      }
#ifdef HYTEG_BUILD_WITH_ADIOS2
      {
         if ( bp4_ )
         {
            if ( useMyrsInsteadOfTimeStepsAsOutputFrequency_ )
            {
               if ( saveNextOutput_ || currentTimeZero )
               {
                  bp4Output_->write( maxLevel_, stepCounter_ );
               }
            }
            else
            {
               if ( writeFrequencyOutput_ > 0 && stepCounter_ % writeFrequencyOutput_ == 0 )
               {
                  bp4Output_->write( maxLevel_, stepCounter_ );
               }
            }
         }
      }
#endif

      // ##########################
      // #### Write Checkpoint ####
      // ##########################
      if ( writeCheckpoints_ )
      {
         if ( useMyrsInsteadOfTimeStepsAsOutputFrequency_ )
         {
            if ( saveNextCheckpoint_ || currentTimeZero )
            {
               createCheckpoint( fileName_ );
            }
         }
         else
         {
            if ( writeFrequencyCheckpoint_ > 0 && stepCounter_ % writeFrequencyCheckpoint_ == 0 )
            {
               createCheckpoint( fileName_ );
            }
         }
      }
   }

   void interpolateDensity( DensityFunctionType_& rho, P2Function< real_t >& T, uint_t level, hyteg::DoFType flag )
   {
      if constexpr ( DensityIsP2_ )
      {
         if constexpr ( TemperatureMassIsP2_ )
         {
            densityModel_->interpolate( rho, { T }, level, flag );
         }
         else
         {
            WALBERLA_ABORT( "Progated temperature is P1 but the density is P2!" );
         }
      }
      else
      {
         densityModel_->interpolate( rho, { T.getVertexDoFFunction() }, level, flag );
      }
   }

   void interpolateViscosity( ViscosityFunctionType_& eta, P2Function< real_t >& T, uint_t level, hyteg::DoFType flag )
   {
      if constexpr ( ViscosityIsP2_ )
      {
         if constexpr ( TemperatureMassIsP2_ )
         {
            viscosityModel_->interpolate( eta, { T }, level, flag );
         }
         else
         {
            WALBERLA_ABORT( "Progated temperature is P1 but the viscosity is P2!" );
         }
      }
      else
      {
         viscosityModel_->interpolate( eta, { T.getVertexDoFFunction() }, level, flag );
      }
   }

   void enforceTemperatureBoundaryConditions( const P2Function< real_t >& T, uint_t level )
   {
      initialTemperatureModel_->interpolate( T, level, dirichletFlag_ );
   }

   void enforceVelocityBoundaryConditions( const P2VectorFunction< real_t >& u, uint_t level )
   {
      if ( CMBType_ == DirichletBoundary )
      {
         setVelocityUIDZero( u, level, "cmb" );
      }

      if ( SurfaceType_ == DirichletBoundary )
      {
         if ( usePlates_ )
         {
            updatePlateVelocities( u, level );
         }
         else
         {
            setVelocityUIDZero( u, level, "surface" );
         }
      }

      projection_->project( u, level, velocityProjectionFlag_ );
   }

   void enforcePressureBoundaryConditions( const P1Function< real_t >& p, uint_t level )
   {
      hyteg::projectPressureMean( p, level );
   }

   void enforceSaddlePointBoundaryConditions( const P2P1TaylorHoodFunction< real_t >& up, uint_t level )
   {
      enforceVelocityBoundaryConditions( up.uvw(), level );
      enforcePressureBoundaryConditions( up.p(), level );
   }

   // IMPORTANT: This assumes we do not need the coarser levels of T for multigrid
   // Otherwise you would have to set the boundary to zero on coarser levels first
   // If we do not do this on P1, we can get overshoots/undershoots in the temperature.
   // This can happen e.g. for steep gradients at plume/slab borders.
   // TODO: Can this be improved?
   void propagateTemperatureToCoarserLevel( P2Function< real_t >& T, uint_t maxLevel, uint_t minLevel )
   {
      std::shared_ptr< P2Function< real_t > > temporaryRHS;

      if ( !lowMemoryMode_ )
      {
         temporaryRHS = tmpTemp_;
      }
      else
      {
         temporaryRHS = getTemporaryFunction< P2Function< real_t > >( storage_, minLevel_, maxLevel_ );

         temporaryRHS->setBoundaryCondition( BC_.bcTemperature_ );
      }

      if constexpr ( TemperatureMassIsP2_ )
      {
         if constexpr ( std::is_same< TemperatureRestrictionOperatorType_, P2toP2InjectionRestriction >::value )
         {
            for ( uint_t level = maxLevel; level > minLevel; level-- )
            {
               temperatureRestriction_->restrict( T, level, allFlag_ );
               enforceTemperatureBoundaryConditions( T, level - 1 );
            }
         }
         else
         {
            temperatureMass_->apply( T, *temporaryRHS, maxLevel, allFlag_, Replace );
            for ( uint_t level = maxLevel; level > minLevel; level-- )
            {
               temperatureRestriction_->restrict( *temporaryRHS, level, allFlag_ );
               enforceTemperatureBoundaryConditions( T, level - 1 );
               temperatureMassSolver_->solve( *temperatureMass_, T, *temporaryRHS, level - 1 );
            }
         }
      }
      else
      {
         temperatureP2ToP1Mass_->apply( T, temporaryRHS->getVertexDoFFunction(), maxLevel, allFlag_, Replace );
         for ( uint_t level = maxLevel; level > minLevel; level-- )
         {
            temperatureRestriction_->restrict( temporaryRHS->getVertexDoFFunction(), level, allFlag_ );
            enforceTemperatureBoundaryConditions( T, level - 1 );
            temperatureMassSolver_->solve(
                *temperatureMass_, T.getVertexDoFFunction(), temporaryRHS->getVertexDoFFunction(), level - 1 );
         }
      }
   }

   void calculateSaddlePointRHS()
   {
      T_d_->interpolate( dynamicTemperature_, { T_->getState( 0 ) }, maxLevel_, allFlag_ );
      saddlePointRHSOperator_->apply(
          up_extra_->uvw(), up_->getState( 0 ).p(), *T_d_, rho_->getState( 0 ), *fg_, maxLevel_, solverFlag_, Replace );
   }

   void calculateAdvectionDiffusionRHS()
   {
      advectionDiffusionRHSOperator_->apply( *T_extra_, *h_, maxLevel_, solverFlag_, Replace );
   }

   void updatePlateVelocities( const P2VectorFunction< real_t >& u, uint_t level )
   {
      //needed for capture below
      uint_t coordIdx;
      real_t plateAge = getCurrentPlateAge();

      //function to return plate velocities, copied and adapted from PlateVelocityDemo.cpp.
      std::function< real_t( const Point3D& ) > Velocity = [this, &coordIdx, &plateAge]( const Point3D& x ) {
         terraneo::vec3D coords{ x[0], x[1], x[2] };
         //get velocity at current plate age (intervals of 1Ma)
         terraneo::vec3D velocity = plateOracle_->getPointVelocity( coords, plateAge );

         return velocity[int_c( coordIdx )] / ND_.uRef_ / plateVelocityScaling_;
      };

      for ( coordIdx = 0; coordIdx < u.getDimension(); coordIdx++ )
      {
         //interpolate current plate velocities at the surface
         BC_.interpolateUIDsWithName( Velocity, level, u[coordIdx], { "surface" }, { coordIdx }, allFlag_ );
      }
      projection_->project( u, level, dirichletFlag_ );
   }

   void setVelocityUIDZero( const P2VectorFunction< real_t >& u, uint_t level, std::string UIDName )
   {
      for ( uint_t coordIdx = 0; coordIdx < u.getDimension(); coordIdx++ )
      {
         BC_.interpolateUIDsWithName( zeroFct_, level, u[coordIdx], { UIDName }, { coordIdx }, allFlag_ );
      }
   }

   // assumes that up and tmpSaddlePoint_ have the same boundary conditions
   hyteg::Point3D
       calculateSaddlePointResidual( P2P1TaylorHoodFunction< real_t >& up, P2P1TaylorHoodFunction< real_t >& fg, uint_t level )
   {
      std::shared_ptr< P2P1TaylorHoodFunction< real_t > > tmpSaddlePointStep;

      if ( !lowMemoryMode_ )
      {
         tmpSaddlePointStep = tmpSaddlePoint_;
      }
      else
      {
         tmpSaddlePointStep = getTemporaryFunction< P2P1TaylorHoodFunction< real_t > >( storage_, minLevel_, maxLevel_ );

         tmpSaddlePointStep->uvw().component( 0 ).setBoundaryCondition( BC_.bcVelocityX_ );
         tmpSaddlePointStep->uvw().component( 1 ).setBoundaryCondition( BC_.bcVelocityY_ );
         if ( storage_->hasGlobalCells() )
         {
            tmpSaddlePointStep->uvw().component( 2 ).setBoundaryCondition( BC_.bcVelocityZ_ );
         }
      }

      saddlePointOperator_->apply( up, *tmpSaddlePointStep, level, solverFlag_ );
      tmpSaddlePointStep->assign( { 1.0, -1.0 }, { *tmpSaddlePointStep, fg }, level, solverFlag_ );

      real_t residual  = std::sqrt( tmpSaddlePointStep->dotGlobal( *tmpSaddlePointStep, level, solverFlag_ ) );
      real_t residualU = std::sqrt( tmpSaddlePointStep->uvw().dotGlobal( tmpSaddlePointStep->uvw(), level, solverFlag_ ) );
      real_t residualP = std::sqrt( tmpSaddlePointStep->p().dotGlobal( tmpSaddlePointStep->p(), level, solverFlag_ ) );

      return { residual, residualU, residualP };
   }

   // assumes that T and tmpTemp_ have the same boundary conditions
   real_t calculateTempResidual( P2Function< real_t >& T, P2Function< real_t >& h, uint_t level )
   {
      std::shared_ptr< P2Function< real_t > > tmpTempStep;

      if ( !lowMemoryMode_ )
      {
         tmpTempStep = tmpTemp_;
      }
      else
      {
         tmpTempStep = getTemporaryFunction< P2Function< real_t > >( storage_, minLevel_, maxLevel_ );

         tmpTempStep->setBoundaryCondition( BC_.bcTemperature_ );
      }

      advectionDiffusionOperator_->apply( T, *tmpTempStep, level, solverFlag_ );
      tmpTempStep->assign( { 1.0, -1.0 }, { *tmpTempStep, h }, level, solverFlag_ );
      return std::sqrt( tmpTempStep->dotGlobal( *tmpTempStep, level, solverFlag_ ) );
   }

   void createFileVectorCheckpoint( const std::string& filename )
   {
      // file name
      std::stringstream dirName;
      dirName << "./vtk/" << filename << "_FileVectorCheckpoint_" << dim_ << "D_" << stepCounter_ << "_Lvl" << minLevel_ << "-"
              << maxLevel_ << "/";

      auto path = std::filesystem::path( dirName.str() );

      if ( !std::filesystem::exists( path ) || !std::filesystem::is_directory( path ) )
      {
         std::filesystem::create_directory( path );
      }

      // save individual functions from function histories
      up_->saveFunctionsViaFileWritingVector( path.string() );
      T_->saveFunctionsViaFileWritingVector( path.string() );
      rho_->saveFunctionsViaFileWritingVector( path.string() );
      eta_->saveFunctionsViaFileWritingVector( path.string() );

      // add function history state data
      std::vector< std::string > stateDataVector;

      up_->addStateDataToStringVector( stateDataVector );
      T_->addStateDataToStringVector( stateDataVector );
      rho_->addStateDataToStringVector( stateDataVector );
      eta_->addStateDataToStringVector( stateDataVector );

      // add step data
      std::stringstream stepCounterStr;
      stepCounterStr << "stepCounter " << stepCounter_ << ";";
      stateDataVector.push_back( stepCounterStr.str() );

      std::stringstream currentTimeStr;
      currentTimeStr << "currentTime " << std::setprecision( std::numeric_limits< real_t >::max_digits10 + 1 ) << currentTime_
                     << ";";
      stateDataVector.push_back( currentTimeStr.str() );

      std::stringstream currentPlateTimeStr;
      currentPlateTimeStr << "currentPlateTime " << std::setprecision( std::numeric_limits< real_t >::max_digits10 + 1 )
                          << currentPlateTime_ << ";";
      stateDataVector.push_back( currentPlateTimeStr.str() );

      // write to state file
      WALBERLA_ROOT_SECTION()
      {
         std::stringstream fName;
         fName << path.string() << "checkpointState.prm";
         std::ofstream file( fName.str(), std::ios::out );

         if ( file.is_open() )
         {
            file << "Parameters\n"
                 << "{\n";
            for ( auto& s : stateDataVector )
            {
               file << "   " << s << "\n";
            }
            file << "}";
         }
         else
         {
            WALBERLA_ABORT( "Error creating file!" );
         }

         file.close();
      }

      // add copy of the parameter file
      WALBERLA_ROOT_SECTION()
      {
         std::stringstream fName;
         fName << path.string() << "ParameterFile.prm";
         std::ofstream file( fName.str(), std::ios::out );

         if ( file.is_open() )
         {
            for ( auto& s : parameterFileVector_ )
            {
               file << s << "\n";
            }
         }
         else
         {
            WALBERLA_ABORT( "Error creating file!" );
         }

         file.close();
      }
   }

   void createCheckpoint( const std::string& filename )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Writing checkpoint..." );

#ifdef HYTEG_BUILD_WITH_ADIOS2
      if ( defaultToUsingAdiosCheckpoints_ )
      {
         std::shared_ptr< hyteg::AdiosCheckpointExporter > exporter = std::make_shared< AdiosCheckpointExporter >( "" );

         std::vector< std::string >               userAttributeNames;
         std::vector< adiosHelpers::adiostype_t > userAttributeValues;

         // register function histories
         up_->registerFunctionsToAdiosCheckpointExporter( exporter );
         T_->registerFunctionsToAdiosCheckpointExporter( exporter );
         rho_->registerFunctionsToAdiosCheckpointExporter( exporter );
         eta_->registerFunctionsToAdiosCheckpointExporter( exporter );

         // add function history state data
         up_->addStateDataToUserAttributes( userAttributeNames, userAttributeValues );
         T_->addStateDataToUserAttributes( userAttributeNames, userAttributeValues );
         rho_->addStateDataToUserAttributes( userAttributeNames, userAttributeValues );
         eta_->addStateDataToUserAttributes( userAttributeNames, userAttributeValues );

         // add step data
         userAttributeNames.push_back( "stepCounter_" );
         userAttributeValues.push_back( stepCounter_ );
         userAttributeNames.push_back( "currentTime_" );
         userAttributeValues.push_back( currentTime_ );
         userAttributeNames.push_back( "currentPlateTime_" );
         userAttributeValues.push_back( currentPlateTime_ );

         // add copy of the parameter file
         WALBERLA_ROOT_SECTION()
         {
            userAttributeNames.push_back( "parameterFileVector_" );
            userAttributeValues.push_back( parameterFileVector_ );
         }

         // file name
         std::stringstream fname;
         fname << filename << "_Checkpoint_" << dim_ << "D_" << stepCounter_ << "_Lvl" << minLevel_ << "-" << maxLevel_;

         exporter->storeCheckpoint( "./vtk", fname.str(), userAttributeNames, userAttributeValues );
      }
      else
      {
         createFileVectorCheckpoint( filename );
      }
#else
      createFileVectorCheckpoint( filename );
#endif

      WALBERLA_LOG_INFO_ON_ROOT( "done." );
   }

   void loadFileVectorCheckpoint( const std::string& filename, uint_t step )
   {
      // file name
      std::stringstream dirName;
      dirName << "./vtk/" << filename << "_FileVectorCheckpoint_" << dim_ << "D_" << step << "_Lvl" << minLevel_ << "-"
              << maxLevel_ << "/";

      auto path = std::filesystem::path( dirName.str() );

      if ( !std::filesystem::exists( path ) || !std::filesystem::is_directory( path ) )
      {
         std::stringstream adiosName;
         adiosName << "./vtk/" << filename << "_Checkpoint_" << dim_ << "D_" << step << "_Lvl" << minLevel_ << "-" << maxLevel_
                   << "/";

         auto adiosPath = std::filesystem::path( adiosName.str() );

         if ( std::filesystem::exists( adiosPath ) && std::filesystem::is_directory( adiosPath ) )
         {
            WALBERLA_ABORT(
                "Checkpoint unavailable. Directory not found. However an ADIOS 2 Checkpoint with a fitting name was found. Make sure you build HyTeG with ADIOS 2 if you want to load ADIOS 2 Checkpoints." );
         }

         WALBERLA_ABORT( "Checkpoint unavailable. Directory not found." );
      }

      // config
      WALBERLA_LOG_INFO_ON_ROOT( "Checkpoint was created with the parameterfile:" );
      WALBERLA_ROOT_SECTION()
      {
         std::stringstream parameterFileName;
         parameterFileName << path.string() << "ParameterFile.prm";

         std::ifstream file( parameterFileName.str(), std::ios::in );

         if ( file.is_open() )
         {
            std::string line;
            while ( std::getline( file, line ) )
            {
               WALBERLA_LOG_INFO_ON_ROOT( line );
            }
         }
         else
         {
            WALBERLA_ABORT( "Error opening file!" );
         }

         file.close();
      }

      // restore function histories
      up_->loadFunctionsViaFileWritingVector( path.string() );
      T_->loadFunctionsViaFileWritingVector( path.string() );
      rho_->loadFunctionsViaFileWritingVector( path.string() );
      eta_->loadFunctionsViaFileWritingVector( path.string() );

      // load state file
      std::stringstream fName;
      fName << path.string() << "checkpointState.prm";

      auto tempCfg = std::make_shared< walberla::config::Config >();
      tempCfg->readParameterFile( fName.str().c_str() );
      walberla::config::Config::BlockHandle stateParams = tempCfg->getOneBlock( "Parameters" );

      // restore function history state data
      up_->loadStateDataFromBlockHandle( stateParams );
      T_->loadStateDataFromBlockHandle( stateParams );
      rho_->loadStateDataFromBlockHandle( stateParams );
      eta_->loadStateDataFromBlockHandle( stateParams );

      // restore step data
      stepCounter_      = stateParams.getParameter< uint_t >( "stepCounter" );
      currentTime_      = stateParams.getParameter< real_t >( "currentTime" );
      currentPlateTime_ = stateParams.getParameter< real_t >( "currentPlateTime" );

      // inverse rho and eta
      if constexpr ( CalculateInvRhoAndEta_ )
      {
         for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
         {
            inv_rho_->assign( { real_c( 1 ) }, { rho_->getState( 0 ) }, level, allFlag_ );
            inv_rho_->invertElementwise( level, allFlag_, false );

            inv_rho_extra_->assign( { real_c( 1 ) }, { *rho_extra_ }, level, allFlag_ );
            inv_rho_extra_->invertElementwise( level, allFlag_, false );

            inv_eta_->assign( { real_c( 1 ) }, { eta_->getState( 0 ) }, level, allFlag_ );
            inv_eta_->invertElementwise( level, allFlag_, false );

            inv_eta_extra_->assign( { real_c( 1 ) }, { *eta_extra_ }, level, allFlag_ );
            inv_eta_extra_->invertElementwise( level, allFlag_, false );
         }
      }
   }

   void loadCheckpoint( const std::string& filename, uint_t step )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Loading checkpoint..." );

#ifdef HYTEG_BUILD_WITH_ADIOS2
      // directory names
      std::stringstream fileVectorName;
      fileVectorName << "./vtk/" << filename << "_FileVectorCheckpoint_" << dim_ << "D_" << step << "_Lvl" << minLevel_ << "-"
                     << maxLevel_ << "/";
      auto fileVectorPath      = std::filesystem::path( fileVectorName.str() );
      bool fileVectorAvailable = std::filesystem::exists( fileVectorPath ) && std::filesystem::is_directory( fileVectorPath );

      std::stringstream adiosName;
      adiosName << "./vtk/" << filename << "_Checkpoint_" << dim_ << "D_" << step << "_Lvl" << minLevel_ << "-" << maxLevel_
                << "/";
      auto adiosPath      = std::filesystem::path( adiosName.str() );
      bool adiosAvailable = std::filesystem::exists( adiosPath ) && std::filesystem::is_directory( adiosPath );

      if ( ( fileVectorAvailable && !defaultToUsingAdiosCheckpoints_ ) || ( !adiosAvailable && fileVectorAvailable ) )
      {
         loadFileVectorCheckpoint( filename, step );
      }
      else
      {
         // file name
         std::stringstream fname;
         fname << filename << "_Checkpoint_" << dim_ << "D_" << step << "_Lvl" << minLevel_ << "-" << maxLevel_;

         auto importer = std::make_shared< AdiosCheckpointImporter >( "./vtk", fname.str(), "" );

         WALBERLA_LOG_INFO_ON_ROOT( "Loading checkpoint " << fname.str() << "..." );

         importer->printCheckpointInfo();

         // config
         WALBERLA_LOG_INFO_ON_ROOT( "Checkpoint was created with the parameterfile:" );
         WALBERLA_ROOT_SECTION()
         {
            std::vector< std::string > parameterFileVectorLoaded_ =
                importer->getUserAttributeValue< std::vector< std::string > >( "parameterFileVector_" );

            for ( auto& s : parameterFileVectorLoaded_ )
            {
               WALBERLA_LOG_INFO_ON_ROOT( s );
            }
         }

         // restore function histories
         up_->loadFunctionsFromAdiosCheckpointImporter( importer );
         T_->loadFunctionsFromAdiosCheckpointImporter( importer );
         rho_->loadFunctionsFromAdiosCheckpointImporter( importer );
         eta_->loadFunctionsFromAdiosCheckpointImporter( importer );

         // restore function history state data
         up_->loadStateDataFromUserAttributes( importer );
         T_->loadStateDataFromUserAttributes( importer );
         rho_->loadStateDataFromUserAttributes( importer );
         eta_->loadStateDataFromUserAttributes( importer );

         // restore step data
         stepCounter_      = importer->getUserAttributeValue< uint_t >( "stepCounter_" );
         currentTime_      = importer->getUserAttributeValue< real_t >( "currentTime_" );
         currentPlateTime_ = importer->getUserAttributeValue< real_t >( "currentPlateTime_" );

         // inverse rho and eta
         if constexpr ( CalculateInvRhoAndEta_ )
         {
            for ( uint_t level = minLevel_; level <= maxLevel_; level++ )
            {
               inv_rho_->assign( { real_c( 1 ) }, { rho_->getState( 0 ) }, level, allFlag_ );
               inv_rho_->invertElementwise( level, allFlag_, false );

               inv_rho_extra_->assign( { real_c( 1 ) }, { *rho_extra_ }, level, allFlag_ );
               inv_rho_extra_->invertElementwise( level, allFlag_, false );

               inv_eta_->assign( { real_c( 1 ) }, { eta_->getState( 0 ) }, level, allFlag_ );
               inv_eta_->invertElementwise( level, allFlag_, false );

               inv_eta_extra_->assign( { real_c( 1 ) }, { *eta_extra_ }, level, allFlag_ );
               inv_eta_extra_->invertElementwise( level, allFlag_, false );
            }
         }
      }
#else
      loadFileVectorCheckpoint( filename, step );
#endif

      WALBERLA_LOG_INFO_ON_ROOT( "done." );
   }

   real_t CFLTimestep( P2VectorFunction< real_t >& u_,
                       uint_t                      level,
                       P2VectorFunction< real_t >& tmpVel,
                       bool                        useGlobalEstimation = false )
   {
      if ( useGlobalEstimation )
      {
         real_t vMax = hyteg::velocityMaxMagnitude( u_, tmpVel.component( 0 ), tmpVel.component( 1 ), level, allFlag_ );

         return CFL_ / vMax * hMin_;
      }
      else
      {
         if ( u_.getStorage()->hasGlobalCells() )
         {
            communication::syncFunctionBetweenPrimitives( u_[0], level );
            communication::syncFunctionBetweenPrimitives( u_[1], level );
            communication::syncFunctionBetweenPrimitives( u_[2], level );

            real_t dt_Local_Min = std::numeric_limits< real_t >::max();

            // loop over all cells
            for ( auto& it : u_.getStorage()->getCells() )
            {
               Cell& cell = *it.second;

               // get hold of the actual data in the velocity
               PrimitiveDataID< FunctionMemory< real_t >, Cell > uXVertexDoFIdx = u_[0].getVertexDoFFunction().getCellDataID();
               PrimitiveDataID< FunctionMemory< real_t >, Cell > uXEdgeDoFIdx   = u_[0].getEdgeDoFFunction().getCellDataID();
               real_t* uXVertexData = cell.getData( uXVertexDoFIdx )->getPointer( level );
               real_t* uXEdgeData   = cell.getData( uXEdgeDoFIdx )->getPointer( level );

               PrimitiveDataID< FunctionMemory< real_t >, Cell > uYVertexDoFIdx = u_[1].getVertexDoFFunction().getCellDataID();
               PrimitiveDataID< FunctionMemory< real_t >, Cell > uYEdgeDoFIdx   = u_[1].getEdgeDoFFunction().getCellDataID();
               real_t* uYVertexData = cell.getData( uYVertexDoFIdx )->getPointer( level );
               real_t* uYEdgeData   = cell.getData( uYEdgeDoFIdx )->getPointer( level );

               PrimitiveDataID< FunctionMemory< real_t >, Cell > uZVertexDoFIdx = u_[2].getVertexDoFFunction().getCellDataID();
               PrimitiveDataID< FunctionMemory< real_t >, Cell > uZEdgeDoFIdx   = u_[2].getEdgeDoFFunction().getCellDataID();
               real_t* uZVertexData = cell.getData( uZVertexDoFIdx )->getPointer( level );
               real_t* uZEdgeData   = cell.getData( uZEdgeDoFIdx )->getPointer( level );

               // loop over micro cells
               for ( const auto& cellType : celldof::allCellTypes )
               {
                  for ( const auto& microCell : celldof::macrocell::Iterator( level, cellType, 0 ) )
                  {
                     // ####################################
                     // #### calculate element diameter ####
                     // ####################################

                     std::array< indexing::Index, 4 > verts =
                         celldof::macrocell::getMicroVerticesFromMicroCell( microCell, cellType );
                     std::array< Point3D, 4 > coords_micro;
                     for ( uint_t k = 0; k < 4; ++k )
                     {
                        coords_micro[k] = vertexdof::macrocell::coordinateFromIndex( level, cell, verts[k] );
                     }

                     Point3D pos0( coords_micro[0][0], coords_micro[0][1], coords_micro[0][2] );
                     Point3D pos1( coords_micro[1][0], coords_micro[1][1], coords_micro[1][2] );
                     Point3D pos2( coords_micro[2][0], coords_micro[2][1], coords_micro[2][2] );
                     Point3D pos3( coords_micro[3][0], coords_micro[3][1], coords_micro[3][2] );

                     // Compute diameter of the as double the circumsphere radius.
                     Point3D ad = pos0 - pos3;
                     Point3D bd = pos1 - pos3;
                     Point3D cd = pos2 - pos3;

                     real_t a       = ad.norm();
                     real_t b       = bd.norm();
                     real_t c       = cd.norm();
                     real_t a_tilde = ( pos1 - pos2 ).norm();
                     real_t b_tilde = ( pos2 - pos0 ).norm();
                     real_t c_tilde = ( pos0 - pos1 ).norm();

                     // heron type formula (from Cayley Menger determinant)
                     real_t s = ( a * a_tilde + b * b_tilde + c * c_tilde ) / real_c( 2 );

                     // calculate volume
                     real_t vol = std::abs( ad.dot( bd.cross( cd ) ) ) / real_c( 6 );

                     // calculate diam
                     real_t diam =
                         real_c( 2 ) * ( std::sqrt( s * ( s - a * a_tilde ) * ( s - b * b_tilde ) * ( s - c * c_tilde ) ) /
                                         ( real_c( 6 ) * vol ) );

                     // ################################################
                     // #### calculate element velocityMaxMagnitude ####
                     // ################################################

                     // obtain data indices of dofs associated with micro-cell
                     std::array< uint_t, 4 > vertexDoFIndices;
                     vertexdof::getVertexDoFDataIndicesFromMicroCell( microCell, cellType, level, vertexDoFIndices );

                     std::array< uint_t, 6 > edgeDoFIndices;
                     edgedof::getEdgeDoFDataIndicesFromMicroCellFEniCSOrdering( microCell, cellType, level, edgeDoFIndices );

                     real_t localVMaxMag;

                     for ( int k = 0; k < 4; k++ )
                     {
                        Point3D velocityVertex( uXVertexData[vertexDoFIndices[uint_c( k )]],
                                                uYVertexData[vertexDoFIndices[uint_c( k )]],
                                                uZVertexData[vertexDoFIndices[uint_c( k )]] );

                        if ( k == 0 )
                        {
                           localVMaxMag = velocityVertex.norm();
                        }
                        else if ( velocityVertex.norm() > localVMaxMag )
                        {
                           localVMaxMag = velocityVertex.norm();
                        }
                     }
                     for ( int k = 4; k < 10; k++ )
                     {
                        Point3D velocityEdge( uXEdgeData[edgeDoFIndices[uint_c( k - 4 )]],
                                              uYEdgeData[edgeDoFIndices[uint_c( k - 4 )]],
                                              uZEdgeData[edgeDoFIndices[uint_c( k - 4 )]] );

                        if ( velocityEdge.norm() > localVMaxMag )
                        {
                           localVMaxMag = velocityEdge.norm();
                        }
                     }

                     real_t dt_Local = CFL_ / localVMaxMag * diam;

                     if ( dt_Local < dt_Local_Min )
                     {
                        dt_Local_Min = dt_Local;
                     }
                  }
               }
            }

            walberla::mpi::allReduceInplace( dt_Local_Min, walberla::mpi::MIN );

            return dt_Local_Min;
         }
         else
         {
            communication::syncFunctionBetweenPrimitives( u_[0], level );
            communication::syncFunctionBetweenPrimitives( u_[1], level );

            real_t dt_Local_Min = std::numeric_limits< real_t >::max();

            // loop over all faces
            for ( auto& it : u_.getStorage()->getFaces() )
            {
               Face& face = *it.second;

               // get hold of the actual data in the velocity
               PrimitiveDataID< FunctionMemory< real_t >, Face > uXVertexDoFIdx = u_[0].getVertexDoFFunction().getFaceDataID();
               PrimitiveDataID< FunctionMemory< real_t >, Face > uXEdgeDoFIdx   = u_[0].getEdgeDoFFunction().getFaceDataID();
               real_t* uXVertexData = face.getData( uXVertexDoFIdx )->getPointer( level );
               real_t* uXEdgeData   = face.getData( uXEdgeDoFIdx )->getPointer( level );

               PrimitiveDataID< FunctionMemory< real_t >, Face > uYVertexDoFIdx = u_[1].getVertexDoFFunction().getFaceDataID();
               PrimitiveDataID< FunctionMemory< real_t >, Face > uYEdgeDoFIdx   = u_[1].getEdgeDoFFunction().getFaceDataID();
               real_t* uYVertexData = face.getData( uYVertexDoFIdx )->getPointer( level );
               real_t* uYEdgeData   = face.getData( uYEdgeDoFIdx )->getPointer( level );

               // loop over micro faces
               for ( const auto& faceType : facedof::allFaceTypes )
               {
                  for ( const auto& microFace : facedof::macroface::Iterator( level, faceType, 0 ) )
                  {
                     // ####################################
                     // #### calculate element diameter ####
                     // ####################################

                     std::array< indexing::Index, 3 > verts =
                         facedof::macroface::getMicroVerticesFromMicroFace( microFace, faceType );

                     std::array< Point3D, 3 > coords_micro;
                     for ( uint_t k = 0; k < 3; k++ )
                     {
                        coords_micro[k] = vertexdof::macroface::coordinateFromIndex( level, face, verts[k] );
                     }

                     Point3D pos0( coords_micro[0][0], coords_micro[0][1], coords_micro[0][2] );
                     Point3D pos1( coords_micro[1][0], coords_micro[1][1], coords_micro[1][2] );
                     Point3D pos2( coords_micro[2][0], coords_micro[2][1], coords_micro[2][2] );

                     // Compute diameter of the as double the circumcircle radius.
                     Point3D ab = pos1 - pos0;
                     Point3D ac = pos2 - pos0;

                     real_t a = ( pos1 - pos2 ).norm();
                     real_t b = ac.norm();
                     real_t c = ab.norm();

                     // Calculate the volume
                     real_t vol = ab.cross( ac ).norm() / real_c( 2 );

                     // calculate diam
                     real_t diam = real_c( 2 ) * ( a * b * c / ( real_c( 4 ) * vol ) );

                     // ################################################
                     // #### calculate element velocityMaxMagnitude ####
                     // ################################################

                     // obtain data indices of dofs associated with micro-face
                     std::array< uint_t, 3 > vertexDoFIndices;
                     vertexdof::getVertexDoFDataIndicesFromMicroFace( microFace, faceType, level, vertexDoFIndices );

                     std::array< uint_t, 3 > edgeDoFIndices;
                     edgedof::getEdgeDoFDataIndicesFromMicroFaceFEniCSOrdering( microFace, faceType, level, edgeDoFIndices );

                     real_t localVMaxMag;

                     for ( int k = 0; k < 3; k++ )
                     {
                        Point3D velocityVertex( uXVertexData[vertexDoFIndices[uint_c( k )]],
                                                uYVertexData[vertexDoFIndices[uint_c( k )]],
                                                real_c( 0 ) );

                        if ( k == 0 )
                        {
                           localVMaxMag = velocityVertex.norm();
                        }
                        else if ( velocityVertex.norm() > localVMaxMag )
                        {
                           localVMaxMag = velocityVertex.norm();
                        }
                     }
                     for ( int k = 3; k < 6; k++ )
                     {
                        Point3D velocityEdge( uXEdgeData[edgeDoFIndices[uint_c( k - 3 )]],
                                              uYEdgeData[edgeDoFIndices[uint_c( k - 3 )]],
                                              real_c( 0 ) );

                        if ( velocityEdge.norm() > localVMaxMag )
                        {
                           localVMaxMag = velocityEdge.norm();
                        }
                     }

                     real_t dt_Local = CFL_ / localVMaxMag * diam;

                     if ( dt_Local < dt_Local_Min )
                     {
                        dt_Local_Min = dt_Local;
                     }
                  }
               }
            }

            walberla::mpi::allReduceInplace( dt_Local_Min, walberla::mpi::MIN );

            return dt_Local_Min;
         }
      }
   }

   // setter / getter functions
   void setOperators( std::shared_ptr< SaddlePointOperatorType_ >           saddlePointOperator,
                      std::shared_ptr< SaddlePointRHSOperatorType_ >        saddlePointRHSOperator,
                      std::shared_ptr< AdvectionDiffusionOperatorType_ >    advectionDiffusionOperator,
                      std::shared_ptr< AdvectionDiffusionRHSOperatorType_ > advectionDiffusionRHSOperator )
   {
      saddlePointOperator_           = saddlePointOperator;
      saddlePointRHSOperator_        = saddlePointRHSOperator;
      advectionDiffusionOperator_    = advectionDiffusionOperator;
      advectionDiffusionRHSOperator_ = advectionDiffusionRHSOperator;
   }

   void setModels( std::shared_ptr< TemperatureModel< real_t > >                   initialTemperatureModel,
                   std::shared_ptr< TemperatureModel< real_t > >                   referenceTemperatureModel,
                   std::shared_ptr< TemperatureDependentDensityModel< real_t > >   densityModel,
                   std::shared_ptr< TemperatureDependentViscosityModel< real_t > > viscosityModel )
   {
      initialTemperatureModel_   = initialTemperatureModel;
      referenceTemperatureModel_ = referenceTemperatureModel;
      densityModel_              = densityModel;
      viscosityModel_            = viscosityModel;
   }

   void setSolvers( std::shared_ptr< SaddlePointSolver< SaddlePointOperatorType_ > >               saddlePointSolver,
                    std::shared_ptr< AdvectionDiffusionSolver< AdvectionDiffusionOperatorType_ > > advectionDiffusionSolver )
   {
      saddlePointSolver_        = saddlePointSolver;
      advectionDiffusionSolver_ = advectionDiffusionSolver;
   }

   void setOperatorUpdater( std::shared_ptr< OperatorUpdater > operatorUpdater )
   {
      operatorUpdater_ = operatorUpdater;
   }

   // returns the current age in Myrs
   real_t getCurrentAge()
   {
      return currentTime_ * ND_.tRef_ / secondsPerMillionYears_;
   }

   // returns the current age in Myrs
   real_t getCurrentPlateAge()
   {
      return currentPlateTime_ * ND_.tRef_ / secondsPerMillionYears_;
   }

   std::shared_ptr< PrimitiveStorage > getStorage()
   {
      return storage_;
   }
   uint_t getMinLevel()
   {
      return minLevel_;
   }
   uint_t getMaxLevel()
   {
      return maxLevel_;
   }
   bool getLowMemoryMode()
   {
      return lowMemoryMode_;
   }
   NondimensionalisationParameters& getNondimensionalisation()
   {
      return ND_;
   }
   walberla::config::Config::BlockHandle& getParameters()
   {
      return parameters_;
   }
   ViscosityFunctionType_& getEta()
   {
      return eta_->getFunctionByIndex( 0 );
   }
   DensityFunctionType_& getRho()
   {
      return rho_->getFunctionByIndex( 0 );
   }
   P2P1TaylorHoodFunction< real_t >& getUp()
   {
      return up_->getFunctionByIndex( 0 );
   }
   P2Function< real_t >& getT()
   {
      return T_->getFunctionByIndex( 0 );
   }
   P2P1TaylorHoodFunction< real_t >& getUpExtra()
   {
      return *up_extra_;
   }
   P2P1TaylorHoodFunction< real_t >& getFg()
   {
      return *fg_;
   }
   P2Function< real_t >& getH()
   {
      return *h_;
   }
   P2Function< real_t >& getTExtra()
   {
      return *T_extra_;
   }
   ViscosityFunctionType_& getEtaExtra()
   {
      return *eta_extra_;
   }
   DensityFunctionType_& getRhoExtra()
   {
      return *rho_extra_;
   }
   std::shared_ptr< ViscosityFunctionType_ > getInvEtaExtra()
   {
      return inv_eta_extra_;
   }
   std::shared_ptr< DensityFunctionType_ > getInvRhoExtra()
   {
      return inv_rho_extra_;
   }
   std::shared_ptr< ViscosityFunctionType_ > getInvEta()
   {
      return inv_eta_;
   }
   std::shared_ptr< DensityFunctionType_ > getInvRho()
   {
      return inv_rho_;
   }
   std::function< bool( const Point3D& x ) >& getSurfaceFct()
   {
      return surfaceFct_;
   }
   std::function< bool( const Point3D& x ) >& getCMBFct()
   {
      return CMBFct_;
   }
   std::shared_ptr< P2ProjectNormalOperator > getVelocityProjection()
   {
      return projection_;
   }
   std::shared_ptr< FunctionHistory< P2Function< real_t >, real_t > >& getTHistory()
   {
      return T_;
   }
   std::shared_ptr< FunctionHistory< P2P1TaylorHoodFunction< real_t >, real_t > >& getUpHistory()
   {
      return up_;
   }
   ConvectionBC& getBC()
   {
      return BC_;
   }
   real_t getSUPG_scaling()
   {
      return SUPG_scaling_;
   }
   real_t getConst_H()
   {
      return const_H_;
   }
   real_t getConst_alpha()
   {
      return const_alpha_;
   }
   real_t getConst_K_T()
   {
      return const_K_T_;
   }
   real_t getConst_k()
   {
      return const_k_;
   }
   real_t getConst_C_p()
   {
      return const_C_p_;
   }
   uint_t getBDFOrder()
   {
      return BDFOrder_;
   }
   std::function< void( const Point3D& in, Point3D& out ) >& getSurfaceNormal()
   {
      return surfaceNormal_;
   }

   std::shared_ptr< VTKOutput > getVtkOutput()
   {
      return vtkOutput_;
   }
#ifdef HYTEG_BUILD_WITH_ADIOS2
   std::shared_ptr< hyteg::AdiosWriter > getBp4Output()
   {
      return bp4Output_;
   }
#endif

   std::ostream& print( std::ostream& os, uint_t offset = 0 ) const
   {
      // clang-format off
      os << std::string( offset, ' ') << "#######################################################"                                                                                                    << "\n";
      os << std::string( offset, ' ') << "############### Mantle Convection Model ###############"                                                                                                    << "\n";
      os << std::string( offset, ' ') << "#######################################################"                                                                                                    << "\n";
      os << std::string( offset, ' ') << "   "    << "------Parameters------"                                                                                                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "prefix_: "                                              << prefix_                                              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "rank_: "                                                << rank_                                                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "size_: "                                                << size_                                                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "blending_: "                                            << blending_                                            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "lowMemoryMode_: "                                       << lowMemoryMode_                                       << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "alwaysDestroyTemporaryFunctions_: "                     << alwaysDestroyTemporaryFunctions_                     << "\n";      
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "vtk_: "                                                 << vtk_                                                 << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "bp4_: "                                                 << bp4_                                                 << "\n"; 
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "writeCheckpoints_: "                                    << writeCheckpoints_                                    << "\n";  
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "defaultToUsingAdiosCheckpoints_: "                      << defaultToUsingAdiosCheckpoints_                      << "\n";  
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "loadCheckpointOnStart_: "                               << loadCheckpointOnStart_                               << "\n"; 
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "loadCheckpointNumber_: "                                << loadCheckpointNumber_                                << "\n";    
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "useMyrsInsteadOfTimeStepsAsOutputFrequency_: "          << useMyrsInsteadOfTimeStepsAsOutputFrequency_          << "\n";   
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "minLevel_: "                                            << minLevel_                                            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "maxLevel_: "                                            << maxLevel_                                            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "dim_: "                                                 << dim_                                                 << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nTan_: "                                                << nTan_                                                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nRad_: "                                                << nRad_                                                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "writeFrequencyOutput_: "                                << writeFrequencyOutput_                                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "writeFrequencyCheckpoint_: "                            << writeFrequencyCheckpoint_                            << "\n";      
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "maximumSaddlePointSolverIterations_: "                  << maximumSaddlePointSolverIterations_                  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nDoFsAdvectionDiffusionGlobal_: "                       << nDoFsAdvectionDiffusionGlobal_                       << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nDoFsSaddlePointGlobal_: "                              << nDoFsSaddlePointGlobal_                              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nDoFsAdvectionDiffusionLocal_: "                        << nDoFsAdvectionDiffusionLocal_                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nDoFsSaddlePointLocal_: "                               << nDoFsSaddlePointLocal_                               << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nDoFsAdvectionDiffusionLocalMax_: "                     << nDoFsAdvectionDiffusionLocalMax_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nDoFsSaddlePointLocalMax_: "                            << nDoFsSaddlePointLocalMax_                            << "\n";            
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nDoFsAdvectionDiffusionLocalAvg_: "                     << nDoFsAdvectionDiffusionLocalAvg_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nDoFsSaddlePointLocalAvg_: "                            << nDoFsSaddlePointLocalAvg_                            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nDoFsAdvectionDiffusionLocalVar_: "                     << nDoFsAdvectionDiffusionLocalVar_                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nDoFsSaddlePointLocalVar_: "                            << nDoFsSaddlePointLocalVar_                            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nMacroElementsGlobal_: "                                << nMacroElementsGlobal_                                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nMacroElementsLocal_: "                                 << nMacroElementsLocal_                                 << "\n";      
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "nMacroElementsLocalMax_: "                              << nMacroElementsLocalMax_                              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "boundaryTolerance_: "                                   << boundaryTolerance_                                   << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "hMin_: "                                                << hMin_                                                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "hMax_: "                                                << hMax_                                                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "absoluteResidualToleranceOuterSaddlePointSolverLoop_: " << absoluteResidualToleranceOuterSaddlePointSolverLoop_ << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "relativeResidualToleranceOuterSaddlePointSolverLoop_: " << relativeResidualToleranceOuterSaddlePointSolverLoop_ << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "randomInitialGuessU_: "                                 << randomInitialGuessU_                                 << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "randomInitialGuessP_: "                                 << randomInitialGuessP_                                 << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "surfaceBoundaryType_: "                                 << surfaceBoundaryType_                                 << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "CMBBoundaryType_: "                                     << CMBBoundaryType_                                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "SUPG_scaling_: "                                        << SUPG_scaling_                                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "const_H_: "                                             << const_H_                                             << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "const_alpha_: "                                         << const_alpha_                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "const_K_T_: "                                           << const_K_T_                                           << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "const_k_: "                                             << const_k_                                             << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "const_C_p_: "                                           << const_C_p_                                           << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "temperatureMassMaxIterations_: "                        << temperatureMassMaxIterations_                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "temperatureMassRelativeTolerance_: "                    << temperatureMassRelativeTolerance_                    << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "temperatureMassAbsoluteTolerance_: "                    << temperatureMassAbsoluteTolerance_                    << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "temperatureMassPrintInfo_: "                            << temperatureMassPrintInfo_                            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "stepCounter_: "                                         << stepCounter_                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "maxSteps_: "                                            << maxSteps_                                            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "dt_: "                                                  << dt_                                                  << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "useGlobalCFL_: "                                        << useGlobalCFL_                                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "minDt_: "                                               << minDt_                                               << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "maxDt_: "                                               << maxDt_                                               << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "minTimestepMyrs_: "                                     << minTimestepMyrs_                                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "maxTimestepMyrs_: "                                     << maxTimestepMyrs_                                     << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "CFL_: "                                                 << CFL_                                                 << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "currentTime_: "                                         << currentTime_                                         << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "secondsPerMillionYears_: "                              << secondsPerMillionYears_                              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "timeDependent_: "                                       << timeDependent_                                       << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "temperatureExtrapolationOrder_: "                       << temperatureExtrapolationOrder_                       << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "velocityExtrapolationOrder_: "                          << velocityExtrapolationOrder_                          << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "BDFOrder_: "                                            << BDFOrder_                                            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "densityDerivativeOrder_: "                              << densityDerivativeOrder_                              << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "viscosityDerivativeOrder_: "                            << viscosityDerivativeOrder_                            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "MMOC_: "                                                << MMOC_                                                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "usePlates_: "                                           << usePlates_                                           << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "loopPlateAge_: "                                        << loopPlateAge_                                        << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "fileTopologies_: "                                      << fileTopologies_                                      << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "fileReconstructions_: "                                 << fileReconstructions_                                 << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "fileName_: "                                            << fileName_                                            << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "plateVelocityScaling_: "                                << plateVelocityScaling_                                << "\n";
      os << std::string( offset, ' ') << "      " << std::setw( 55 ) << std::left << "currentPlateTime_: "                                    << currentPlateTime_                                    << "\n";
      // clang-format on

      ND_.print( os, offset + 3 );
      os << "\n";
      initialTemperatureModel_->print( os, offset + 3 );
      os << "\n";
      referenceTemperatureModel_->print( os, offset + 3 );
      os << "\n";
      densityModel_->print( os, offset + 3 );
      os << "\n";
      viscosityModel_->print( os, offset + 3 );
      os << "\n";
      saddlePointSolver_->print( os, offset + 3 );
      os << "\n";
      advectionDiffusionSolver_->print( os, offset + 3 );

      return os;
   }

   ~MantleConvectionModel(){};

 protected:
   ///////////////////////////////////////////////////////////////////
   //////////////////// general & mesh parameters ////////////////////
   ///////////////////////////////////////////////////////////////////

   const uint_t rank_;                  // MPI rank of the process
   const uint_t size_;                  // total number of MPI processes
   std::string  prefix_;                // prefix used for ND_ and general parameters
   const bool   blending_;              // true if blending is used
   bool         vtk_;                   // output simulation to vtk file
   bool         bp4_;                   // output simulation to bp4 file, requires Adios2
   bool         writeCheckpoints_;      // true if adios2 checkpoints are written
   bool         loadCheckpointOnStart_; // true if we try to load a checkpoint on startup
   bool
        defaultToUsingAdiosCheckpoints_; // set to false if you want to save or preferably load FileVector checkpoints even if adios2 is available
   bool lowMemoryMode_;                  // low memory mode tries to keep the required memory as low as possible (but is slower)
   bool
        useMyrsInsteadOfTimeStepsAsOutputFrequency_; // true of checkpoints and vtk/bp4 files should be created dependent on the simulation time in Myr instead of the number of timesteps
   bool saveNextCheckpoint_;                         // state variable used if useMyrsInsteadOfTimeStepsAsOutputFrequency_ = true
   bool saveNextOutput_;                             // state variable used if useMyrsInsteadOfTimeStepsAsOutputFrequency_ = true
   bool
                                         alwaysDestroyTemporaryFunctions_; // true if tm functions generated via the TempFunctionManager (lowMemoryMode_ = true) should always be deleted after use
   bool                                  checkpointNotLoaded_; // true if we have not tried to load a checkpoint at the startup
   walberla::config::Config::BlockHandle parameters_;          // parameter block
   std::vector< std::string >            parameterFileVector_; // line by line representation of the parameter file

   uint_t minLevel_;                           // minimum grid refinement level
   uint_t maxLevel_;                           // maximum grid refinement level
   uint_t dim_;                                // dimension of the model ( 2 or 3 are supported )
   uint_t nTan_;                               // number of mesh nodes along spherical diamond edge
   uint_t nRad_;                               // number of radial layers in the mesh
   uint_t writeFrequencyOutput_;               // vtk and bp4 output write frequency
   uint_t writeFrequencyCheckpoint_;           // checkpoint output write frequency
   uint_t loadCheckpointNumber_;               // number of the step we try to load via checkpoint
   uint_t maximumSaddlePointSolverIterations_; // maximum number of saddle point solver iterations (outer loop)
   uint_t nDoFsAdvectionDiffusionGlobal_;      // number of global DoFs of a temperature function on maxLevel_
   uint_t nDoFsSaddlePointGlobal_;             // number of global DoFs of a velocity&pressure composite function on maxLevel_
   uint_t nDoFsAdvectionDiffusionLocal_;       // number of local DoFs of a temperature function on maxLevel_
   uint_t nDoFsSaddlePointLocal_;              // number of local DoFs of a velocity&pressure composite function on maxLevel_
   uint_t nDoFsAdvectionDiffusionLocalMax_;    // maximum number of local DoFs of a temperature function on maxLevel_
   uint_t nDoFsSaddlePointLocalMax_;        // maximum number of local DoFs of a velocity&pressure composite function on maxLevel_
   uint_t nMacroElementsGlobal_;            // number of global macro elements
   uint_t nMacroElementsLocal_;             // number of local macro elements
   uint_t nMacroElementsLocalMax_;          // maximum number of local macro elements
   real_t nDoFsAdvectionDiffusionLocalAvg_; // average number of local DoFs of a temperature function on maxLevel_
   real_t nDoFsSaddlePointLocalAvg_;        // average number of local DoFs of a velocity&pressure composite function on maxLevel_
   real_t nDoFsAdvectionDiffusionLocalVar_; // variance sigma^2 of the number of local DoFs of a temperature function on maxLevel_
   real_t
       nDoFsSaddlePointLocalVar_; // variance sigma^2 of the number of local DoFs of a velocity&pressure composite function on maxLevel_
   real_t boundaryTolerance_; // tolerance for setting the boundary flags and evaluating the surface and CMB functions
   real_t hMin_;              // minimum edge length in the mesh
   real_t hMax_;              // maximum edge length in the mesh
   real_t absoluteResidualToleranceOuterSaddlePointSolverLoop_; // absolute tolerance of the outer saddle point solver loop
   real_t relativeResidualToleranceOuterSaddlePointSolverLoop_; // relative tolerance of the outer saddle point solver loop
   bool   randomInitialGuessU_;            // true if we initialise U with random values otherwise a zero initial guess is used
   bool   randomInitialGuessP_;            // true if we initialise P with random values otherwise a zero initial guess is used
   hyteg::DoFType SurfaceType_;            // surface DoFType
   hyteg::DoFType CMBType_;                // CMB DofType
   hyteg::DoFType solverFlag_;             // solver flag
   hyteg::DoFType dirichletFlag_;          // dirichlet flag
   hyteg::DoFType allFlag_;                // all flag
   hyteg::DoFType velocityProjectionFlag_; // velocity projection flag

   real_t SUPG_scaling_; // variable scaling for the SUPG terms (default 1.0)

   // values that we set to a constant in our model
   real_t const_H_;     // nondimensional internal heating             , set to 1.0 for H_Ref
   real_t const_alpha_; // nondimensional thermal expansivity          , set to 1.0 for alpha_Ref
   real_t const_K_T_;   // nondimensional isothermal bulk modulus value, set to 1.0 for K_T_Ref
   real_t const_k_;     // nondimensional thermal conductivity         , set to 1.0 for k_Ref
   real_t const_C_p_;   // nondimensional specific heat capacity       , set to 1.0 for C_p__Ref

   uint_t temperatureMassMaxIterations_; // maximum number of iterations of the temperature propagation solver to coarser levels
   real_t temperatureMassRelativeTolerance_; // relative tolerance of the temperature mass solver
   real_t temperatureMassAbsoluteTolerance_; // absolute tolerance the temperature mass solver
   bool   temperatureMassPrintInfo_;         // print temperature mass solver output

   uint_t stepCounter_;            // current simulation time step
   uint_t maxSteps_;               // maximum number of time steps
   real_t dt_;                     // current time step
   real_t minDt_;                  // minimum value for dt_, dep. on minTimestepMyrs_
   real_t maxDt_;                  // maximum value for dt_, dep. on maxTimestepMyrs_
   real_t minTimestepMyrs_;        // minimum time step size in Myrs
   real_t maxTimestepMyrs_;        // maximum time step size in Myrs
   real_t CFL_;                    // CFL Constant
   real_t currentTime_;            // nondimensional current time
   real_t secondsPerMillionYears_; // seconds in a million years
   bool
                         useGlobalCFL_; // true if we use a global CFL estimation to determine the next time step instead of taking all individual element conditions into account
   static constexpr bool timeDependent_ =
       !( std::is_same< typename AdvectionDiffusionOperatorType_::MassOperatorTypeInternal, hyteg::NoOperator >::value &&
          std::is_same< typename AdvectionDiffusionOperatorType_::MassStabilisationOperatorTypeInternal,
                        hyteg::NoOperator >::value ); // true if we solve a time dependent system

   ConvectionBC BC_;                  // boundary conditions and uids
   uint_t       surfaceBoundaryType_; // surface boundary type, 1 = Dirichlet, 2 = Neumann, 3 = FreeSlip
   uint_t       CMBBoundaryType_;     // cmb boundary type, 1 = Dirichlet, 2 = Neumann, 3 = FreeSlip

   uint_t temperatureExtrapolationOrder_; // order of the temperature extrapolation
   uint_t velocityExtrapolationOrder_;    // order of the velocity extrapolation
   uint_t BDFOrder_;                      // order of the BDF Scheme
   uint_t densityDerivativeOrder_;        // order of the density derivative approximation (needed for projected density approx)
   uint_t viscosityDerivativeOrder_;      // order of the viscosity derivative approximation (currently not used)

   static constexpr bool MMOC_ = std::is_same< typename AdvectionDiffusionOperatorType_::AdvectionOperatorTypeInternal,
                                               hyteg::NoOperator >::value; // true if we use the MMOC method

   static constexpr bool TemperatureMassIsP2_ =
       std::is_same< typename TemperatureMassOperatorType_::dstType,
                     hyteg::P2Function< real_t > >::value; // true if the temperature is propagated two lower levels via a P2Mass
   static constexpr bool DensityIsP2_ =
       std::is_same< DensityFunctionType_,
                     hyteg::P2Function< real_t > >::value; // true if the denisty is represented via a P2 FEM function
   static constexpr bool ViscosityIsP2_ =
       std::is_same< ViscosityFunctionType_,
                     hyteg::P2Function< real_t > >::value; // true if the viscosity is represented via a P2 FEM function

   bool        usePlates_;      // true if we use the plate velocities from the json file to update the surface boundary
   bool        loopPlateAge_;   // true if we loop the plate boundaries (in a modulo fashion) after the max age saved in the file
   std::string fileTopologies_; // filepath for the topologies file
   std::string fileReconstructions_;  // filepath for the plate reconstruction file
   std::string fileName_;             // filename for the output of files
   real_t      plateVelocityScaling_; // plate velocity scaling
   real_t      currentPlateTime_;     // nondimensional current plate

   bool loadStorageFile_; // true if we try to load the setupPrimitiveStorage from fileSetupPrimitiveStorage_ on startup
   std::string
       fileSetupPrimitiveStorage_; // filepath for the setup primitive storage file to load, only used when loadStorageFile_ = true

   ///////////////////////////////////////////////////////
   //////////////////// FEM Functions ////////////////////
   ///////////////////////////////////////////////////////

   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > up_extra_;       // extrapolation of up_ to the next time step
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > fg_;             // right hand side of the saddle point system
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > tmpSaddlePoint_; // temporary saddle point function

   std::shared_ptr< P2Function< real_t > >   T_d_;               // current dynamic temperature
   std::shared_ptr< P2Function< real_t > >   T_extra_;           // extrapolation of T_ to the next time step
   std::shared_ptr< P2Function< real_t > >   h_;                 // right hand side of the temperature system
   std::shared_ptr< P2Function< real_t > >   tmpTemp_;           // temporary temperature function
   std::shared_ptr< P2Function< real_t > >   tmpTemp2_;          // second temporary temperature function
   std::shared_ptr< P2Function< real_t > >   mmocDummyFunction_; // dummy function for MMOC second order
   std::shared_ptr< ViscosityFunctionType_ > eta_extra_;         // viscosity extrapolation
   std::shared_ptr< DensityFunctionType_ >   rho_extra_;         // density extrapolation

   std::shared_ptr< ViscosityFunctionType_ > inv_eta_;       // inv viscosity
   std::shared_ptr< DensityFunctionType_ >   inv_rho_;       // inv density
   std::shared_ptr< ViscosityFunctionType_ > inv_eta_extra_; // inv viscosity extrapolation
   std::shared_ptr< DensityFunctionType_ >   inv_rho_extra_; // inv density extrapolation

   std::shared_ptr< FunctionHistory< P2P1TaylorHoodFunction< real_t >, real_t > > up_;  // velocity & pressure function history
   std::shared_ptr< FunctionHistory< P2Function< real_t >, real_t > >             T_;   // temperature function history
   std::shared_ptr< FunctionHistory< DensityFunctionType_, real_t > >             rho_; // density function history
   std::shared_ptr< FunctionHistory< ViscosityFunctionType_, real_t > >           eta_; // viscosity function history

   ///////////////////////////////////////////////////////////////////
   //////////////////// internal member variables ////////////////////
   ///////////////////////////////////////////////////////////////////

   std::shared_ptr< PrimitiveStorage > storage_;   // primitive storage
   NondimensionalisationParameters     ND_;        // Nondimensionalisation parameters
   std::shared_ptr< VTKOutput >        vtkOutput_; // vtk output class for the simulation
#ifdef HYTEG_BUILD_WITH_ADIOS2
   std::shared_ptr< hyteg::AdiosWriter > bp4Output_; // bp4 output class for the simulation
#endif
   std::shared_ptr< walberla::config::Config > cfg_; // walberla config

   std::function< bool( const Point3D& x ) > surfaceFct_; // returns true if x is on the surface (up to boundaryTolerance_)
   std::function< bool( const Point3D& x ) > CMBFct_;     // returns true if x is on the cmb (up to boundaryTolerance_)
   std::function< void( const Point3D& in, Point3D& out ) > surfaceNormal_; // returns the surface and cmb normal direction at in
   std::function< real_t( const Point3D& ) >                zeroFct_;       // return zero
   std::function< real_t( const Point3D& ) >                randFunc_;      // returns a random number between -1 and 1
   std::function< real_t( const hyteg::Point3D& ) >         referenceTemperature_; // evaluates the reference temperature
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) >
       dynamicTemperature_; // evaluates the dynamic temperature T_d = T - T_ref

   std::shared_ptr< P2ProjectNormalOperator > projection_; // free slip normal projection operator

   std::shared_ptr< MMOCTransport< P2Function< real_t > > > mmocTransport_; // mmoc transport

   std::shared_ptr< terraneo::plates::PlateVelocityProvider > plateOracle_; // terraneo plate velocity provider / oracle

   std::shared_ptr< SaddlePointOperatorType_ >           saddlePointOperator_;
   std::shared_ptr< SaddlePointRHSOperatorType_ >        saddlePointRHSOperator_;
   std::shared_ptr< AdvectionDiffusionOperatorType_ >    advectionDiffusionOperator_;
   std::shared_ptr< AdvectionDiffusionRHSOperatorType_ > advectionDiffusionRHSOperator_;

   std::shared_ptr< TemperatureModel< real_t > >                   initialTemperatureModel_;
   std::shared_ptr< TemperatureModel< real_t > >                   referenceTemperatureModel_;
   std::shared_ptr< TemperatureDependentDensityModel< real_t > >   densityModel_;
   std::shared_ptr< TemperatureDependentViscosityModel< real_t > > viscosityModel_;

   std::shared_ptr< SaddlePointSolver< SaddlePointOperatorType_ > >               saddlePointSolver_;
   std::shared_ptr< AdvectionDiffusionSolver< AdvectionDiffusionOperatorType_ > > advectionDiffusionSolver_;

   std::shared_ptr< P2ToP1ElementwiseBlendingKMassOperator >
       temperatureP2ToP1Mass_; // mass operator used to project the temperature down in levels
   std::shared_ptr< TemperatureMassOperatorType_ >
       temperatureMass_; // mass operator used to project the temperature down in levels
   std::shared_ptr< TemperatureRestrictionOperatorType_ >
       temperatureRestriction_; // restriction operator used to project the temperature down in levels
   std::shared_ptr< hyteg::CGSolver< TemperatureMassOperatorType_ > >
       temperatureMassSolver_; // solver used to project the temperature down in levels

   std::shared_ptr< OperatorUpdater >
       operatorUpdater_; // operator updater that can be customised to e.g. recalculate inverseDiagonals when necessary
};

template < class SaddlePointOperatorType_,
           class SaddlePointRHSOperatorType_,
           class AdvectionDiffusionOperatorType_,
           class AdvectionDiffusionRHSOperatorType_,
           class TemperatureMassOperatorType_,
           class TemperatureRestrictionOperatorType_ >
inline std::ostream& operator<<( std::ostream&                                                       os,
                                 const MantleConvectionModel< SaddlePointOperatorType_,
                                                              SaddlePointRHSOperatorType_,
                                                              AdvectionDiffusionOperatorType_,
                                                              AdvectionDiffusionRHSOperatorType_,
                                                              TemperatureMassOperatorType_,
                                                              TemperatureRestrictionOperatorType_ >& mcm )
{
   return mcm.print( os );
}

} // namespace MantleConvection
