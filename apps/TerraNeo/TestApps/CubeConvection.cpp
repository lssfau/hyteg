/*
 * Copyright (c) 2025 Ponsuganth Ilangovan
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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointExporter.hpp"
#include "hyteg/checkpointrestore/ADIOS2/AdiosCheckpointImporter.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/composites/UnsteadyDiffusion.hpp"
#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/forms/P2LinearCombinationForm.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_mass_blending_q4.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP2Conversion.hpp"
#include "hyteg/gridtransferoperators/P2toP1Conversion.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticVectorRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/FGMRESSolver.hpp"
#include "hyteg/solvers/GMRESSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockPreconditioners.hpp"
#include "hyteg_operators/operators/k_divdiv/P2VectorElementwiseKDivdiv.hpp"
#include "hyteg_operators/operators/k_mass/P1ElementwiseKMass.hpp"
#include "hyteg_operators/operators/k_mass/P2ToP1ElementwiseKMass.hpp"
#include "hyteg_operators/operators/mass/P1ElementwiseMass.hpp"
#include "hyteg_operators/operators/mass/P2ElementwiseMass.hpp"
#include "hyteg_operators/operators/terraneo/P2VectorToP1ElementwiseFrozenVelocity.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesConstantOperator.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesEpsilonOperator.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

// #include "hyteg/operatorgeneration/generated/EvaluateViscosityViscoplastic/P2EvaluateViscosityViscoplastic.hpp"
// #include "hyteg/operatorgeneration/generated/FullStokesViscoplastic/P2VectorElementwiseFullStokesViscoplastic.hpp"

#include "coupling_hyteg_convection_particles/MMOCTransport.hpp"
#include "mixed_operator/VectorMassOperator.hpp"
#include "terraneo/helpers/RadialProfiles.hpp"
#include "terraneo/operators/TransportOperatorStd.hpp"
#include "terraneo/utils/NusseltNumberOperator.hpp"

using walberla::real_t;
using walberla::uint_t;

using namespace hyteg;
using namespace terraneo;

using P2ToP1ElementwiseKMass = operatorgeneration::P2ToP1ElementwiseKMass;
using FrozenVelocityOperator = operatorgeneration::P2VectorToP1ElementwiseFrozenVelocity;

namespace hyteg {

namespace solvertemplates {

template < typename StokesOperatorType, typename StokesABlockType, typename StokesSchurOperatorType >
inline std::shared_ptr< Solver< StokesOperatorType > > fgmresMGSolver( const std::shared_ptr< PrimitiveStorage >& storage,
                                                                       uint_t                                     minLevel,
                                                                       uint_t                                     maxLevel,
                                                                       const StokesABlockType&        stokesABlockOperator,
                                                                       const StokesSchurOperatorType& schurOperator )
{
   uint_t ABlockPreSmooth  = 2u;
   uint_t ABlockPostSmooth = 2u;

   auto ABlockSmoother = std::make_shared< ChebyshevSmoother< StokesABlockType > >( storage, minLevel, maxLevel );

   auto uTmp = std::make_shared< typename StokesABlockType::srcType >( "uTmpSolverTemplate", storage, minLevel, maxLevel );
   auto uSpecTmp =
       std::make_shared< typename StokesABlockType::srcType >( "uSpecTmpSolverTemplate", storage, minLevel, maxLevel );

   std::function< real_t( const Point3D& ) > randFuncA = []( const Point3D& ) {
      return walberla::math::realRandom( real_c( -1 ), real_c( 1 ) );
   };

   uTmp->interpolate( { randFuncA, randFuncA, randFuncA }, maxLevel, All );
   uSpecTmp->interpolate( { randFuncA, randFuncA, randFuncA }, maxLevel, All );

   real_t spectralRadius = chebyshev::estimateRadius( stokesABlockOperator, maxLevel, 25u, storage, *uTmp, *uSpecTmp );

   WALBERLA_LOG_INFO_ON_ROOT( "spectralRadius = " << spectralRadius );

   ABlockSmoother->setupCoefficients( 1u, spectralRadius );

   auto ABlockProlongationOperator = std::make_shared< P2toP2QuadraticVectorProlongation >();
   auto ABlockRestrictionOperator  = std::make_shared< P2toP2QuadraticVectorRestriction >();

   std::shared_ptr< Solver< StokesABlockType > > ABlockCoarseGridSolver;

#if defined( HYTEG_BUILD_WITH_PETSC )
   ABlockCoarseGridSolver = std::make_shared< PETScLUSolver< StokesABlockType > >( storage, minLevel );
#else
   ABlockCoarseGridSolver = std::make_shared< MinResSolver< StokesABlockType > >( storage, minLevel, minLevel, 100u );
#endif

   auto ABlockMultigridSolver = std::make_shared< GeometricMultigridSolver< StokesABlockType > >( storage,
                                                                                                  ABlockSmoother,
                                                                                                  ABlockCoarseGridSolver,
                                                                                                  ABlockRestrictionOperator,
                                                                                                  ABlockProlongationOperator,
                                                                                                  minLevel,
                                                                                                  maxLevel,
                                                                                                  ABlockPreSmooth,
                                                                                                  ABlockPostSmooth,
                                                                                                  0,
                                                                                                  CycleType::VCYCLE );

   auto identityPreconditioner = std::make_shared< IdentityPreconditioner< StokesABlockType > >();

   uint_t SchurOuterIter = 500u;
   real_t SchurOuterTol  = 1e-12;

   auto SchurSolver =
       std::make_shared< CGSolver< StokesSchurOperatorType > >( storage, minLevel, maxLevel, SchurOuterIter, SchurOuterTol );

   auto blockPreconditioner =
       std::make_shared< BlockFactorisationPreconditioner< StokesOperatorType, StokesABlockType, StokesSchurOperatorType > >(
           storage, minLevel, maxLevel, schurOperator, ABlockMultigridSolver, SchurSolver, 1u );

   uint_t fGMRESOuterIter = 10u;
   real_t fGMRESTol       = 1e-6;

   auto finalStokesSolver = std::make_shared< FGMRESSolver< StokesOperatorType > >(
       storage, minLevel, maxLevel, fGMRESOuterIter, fGMRESTol, fGMRESTol, blockPreconditioner, 50u, fGMRESTol );
   finalStokesSolver->setPrintInfo( true );

   return finalStokesSolver;
}

} // namespace solvertemplates

// namespace operatorgeneration {

// using P2P1StokesViscoplastic =
//     detail::P2P1StokesNonlinViscOperatorTemplate< operatorgeneration::P2VectorElementwiseFullStokesViscoplastic,
//                                                   operatorgeneration::P1ToP2GradientOperator,
//                                                   operatorgeneration::P2ToP1DivergenceOperator >;
// } // namespace operatorgeneration

struct ParameterContainer
{
   bool verbose = true;

   real_t rMin = 1.22, rMax = 2.22;

   uint_t maxTimeSteps = 1000, vtkWriteFrequency = 1U;

   bool MMOC = true, SUPG = false, compressible = true, shearHeating = true;

   real_t Ra = 1e5, Di = 0.5, T0 = 0.091, diffusivity = 1.0, cflMax = 0.75, AiniPerturb = 0.1;

   real_t rho0 = 1.0, alpha = 1.0, cpr = 1.0, cvr = 1.0, grueneisen = 1.0, alphabar = 1.0, cpbar = 1.0, chibar = 1.0, k_ = 1.0;

   real_t minresRelTol = 1e-4, minresAbsTol = 1e-8, gmresTol = 1e-5;
   uint_t minresIter = 1000U, gmresIter = 1000U;

   uint_t nsCalcFreq = 10U;
};

const real_t boundaryMarkerThreshold = 1e-6;

const real_t hLmaxCuboid__ = 1.0;
const real_t vLmaxCuboid__ = 1.0;

std::function< bool( const Point3D& ) > bottomMarker = []( const Point3D& x ) {
   if ( std::abs( x[2] ) < boundaryMarkerThreshold )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > rightMarker = []( const Point3D& x ) {
   if ( std::abs( x[0] - hLmaxCuboid__ ) < boundaryMarkerThreshold )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > leftMarker = []( const Point3D& x ) {
   if ( std::abs( x[0] ) < boundaryMarkerThreshold )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > frontMarker = []( const Point3D& x ) {
   if ( std::abs( x[1] ) < boundaryMarkerThreshold )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > backMarker = []( const Point3D& x ) {
   if ( std::abs( x[1] - hLmaxCuboid__ ) < boundaryMarkerThreshold )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > topMarker = []( const Point3D& x ) {
   if ( std::abs( x[2] - vLmaxCuboid__ ) < boundaryMarkerThreshold )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > BRMarker = []( const Point3D& x ) {
   if ( ( bottomMarker( x ) && rightMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > BLMarker = []( const Point3D& x ) {
   if ( ( bottomMarker( x ) && leftMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > BFMarker = []( const Point3D& x ) {
   if ( ( bottomMarker( x ) && frontMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > BBMarker = []( const Point3D& x ) {
   if ( ( bottomMarker( x ) && backMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > TRMarker = []( const Point3D& x ) {
   if ( ( topMarker( x ) && rightMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > TLMarker = []( const Point3D& x ) {
   if ( ( topMarker( x ) && leftMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > TFMarker = []( const Point3D& x ) {
   if ( ( topMarker( x ) && frontMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > TBMarker = []( const Point3D& x ) {
   if ( ( topMarker( x ) && backMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > FRMarker = []( const Point3D& x ) {
   if ( ( frontMarker( x ) && rightMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > FLMarker = []( const Point3D& x ) {
   if ( ( frontMarker( x ) && leftMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > BaRMarker = []( const Point3D& x ) {
   if ( ( backMarker( x ) && rightMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > BaLMarker = []( const Point3D& x ) {
   if ( ( backMarker( x ) && leftMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

std::function< bool( const Point3D& ) > cornersMarker = []( const Point3D& x ) {
   if ( ( topMarker( x ) && rightMarker( x ) && frontMarker( x ) ) ||
        ( bottomMarker( x ) && rightMarker( x ) && frontMarker( x ) ) ||
        ( topMarker( x ) && leftMarker( x ) && frontMarker( x ) ) ||
        ( bottomMarker( x ) && leftMarker( x ) && frontMarker( x ) ) ||
        ( topMarker( x ) && rightMarker( x ) && backMarker( x ) ) ||
        ( bottomMarker( x ) && rightMarker( x ) && backMarker( x ) ) ||
        ( topMarker( x ) && leftMarker( x ) && backMarker( x ) ) || ( bottomMarker( x ) && leftMarker( x ) && backMarker( x ) ) )
   {
      return true;
   }
   else
   {
      return false;
   }
};

// typedef operatorgeneration::P2P1StokesFullP0ViscosityOperator StokesOperator;
// // typedef operatorgeneration::P2P1StokesFullP1ViscosityOperator StokesOperator;
// typedef P2P1ElementwiseBlendingStokesOperator      StokesOperatorLinear;
// typedef operatorgeneration::P2P1StokesViscoplastic StokesOperatorNonlinear;

// typedef StrongFreeSlipWrapper< StokesOperator, P2ProjectNormalOperator, true > StokesOperatorFS;

using StokesOperator_T = operatorgeneration::P2P1StokesEpsilonOperator;
using SchurOperator_T  = operatorgeneration::P1ElementwiseKMass;

enum BoundaryMarkers
{
   Bottom = 23,
   Right,
   Left,
   Top,
   Front,
   Back,
   BRCorner,
   BLCorner,
   BFCorner,
   BBCorner,
   TRCorner,
   TLCorner,
   TFCorner,
   TBCorner,
   FRCorner,
   FLCorner,
   BaRCorner,
   BaLCorner,
   PointCorners,
};

real_t viscosityFunction( const Point3D& x, const std::vector< real_t >& Temperature )
{
   auto   radius = x[2];
   real_t retVal = 1;
   real_t rMax   = 2.22;

   real_t activationEnergy     = real_c( 3 );
   real_t depthViscosityFactor = real_c( 1.5 );
   real_t viscosityLowerBound  = real_c( 1e20 );
   real_t viscosityUpperBound  = real_c( 1e24 );
   real_t c1                   = real_c( 0.05 );
   real_t c2                   = real_c( 2.05 );

   //set background viscosity as radial profile or constant value, and non-dimensionalise
   // if ( haveViscosityProfile )
   // {
   //    retVal = linearInterpolateBetween( viscosityProfile, radius );
   // }
   // else
   // {
   //    retVal = parameters.referenceViscosity;
   // }

   //scale background viscosity by temperature- and depth-dependent factors
   //depth-dependent factor counteracts the decrease in viscosity due to increasing temperature with depth
   // if ( tempDependentViscosity )
   // {
   // switch ( tempDependentViscosityType )
   // {
   // //Frank–Kamenetskii type 1
   // case 0: {
   retVal *= std::exp( -activationEnergy * ( Temperature[0] ) + depthViscosityFactor * ( rMax - radius ) / ( rMax ) );
   //    break;
   // }
   // //Frank–Kamenetskii type 2
   // case 1: {
   //    retVal *= std::exp( activationEnergy * ( real_c( 0.5 ) - Temperature[0] ) +
   //                        depthViscosityFactor * ( parameters.rMax - radius ) / ( parameters.rMax ) );
   //    break;
   // }
   // //Arrhenius type
   // case 2: {

   // WALBERLA_LOG_INFO_ON_ROOT( "temp-dep viscosity : Arrhenius" );
   // WALBERLA_LOG_INFO_ON_ROOT( "Temperature : " );
   // WALBERLA_LOG_INFO_ON_ROOT( Temperature[0] );

   // retVal *= std::exp( activationEnergy * ( ( real_c( 1 ) / ( Temperature[0] + c1 ) ) - c2 ) +
   //                           depthViscosityFactor * ( rMax - radius ) / ( rMax ) );

   //       break;
   //    }
   //    //Frank–Kamenetskii type 1
   //    default: {
   //       retVal *= std::exp( -activationEnergy * ( Temperature[0] ) +
   //                           depthViscosityFactor * ( parameters.rMax - radius ) / ( parameters.rMax ) );
   //       break;
   //    }
   //    }
   //    // WALBERLA_LOG_INFO_ON_ROOT( " viscosity before non - D : " << retVal );
   //    //  impose min viscosity
   if ( retVal < 1 )
   {
      retVal = 1;
   }

   //    //impose max viscosity
   if ( retVal > 100 )
   {
      retVal = 100;
   }
   // // }
   // if ( weakZone )
   // {
   //    // for the sake of testing out, let's define 3 random points as weak zones
   //    // gotta test with and without plate velocities

   //    real_t X_1 = 0;
   //    real_t Y_1 = parameters.rMax;

   //    real_t X_2 = parameters.rMax * std::cos( ( walberla::math::pi ) / 6 );
   //    real_t Y_2 = ( -1 ) * parameters.rMax * std::sin( ( walberla::math::pi ) / 6 );

   //    real_t X_3 = ( -1 ) * parameters.rMax * std::cos( ( walberla::math::pi ) / 4 );
   //    real_t Y_3 = ( -1 ) * parameters.rMax * std::sin( ( walberla::math::pi ) / 4 );

   //    // now let's make sure they have lower viscosity

   //    const real_t rS = 200000 / ( parameters.rSurface );
   //    const real_t rL = 100000 / ( parameters.rSurface );

   //    // looking bad
   //    if ( x[0] > X_1 && x[0] < rS )
   //    {
   //       if ( x[1] > ( Y_1 - rL ) )
   //       {
   //          retVal *= 0.5;
   //       }
   //    }
   //    const real_t rL_2 = 100000 / ( parameters.rSurface );
   //    // looking good
   //    if ( x[0] < X_3 )
   //    {
   //       if ( x[1] < ( Y_3 + rL_2 ) )
   //       {
   //          retVal *= 0.1;
   //       }
   //    }
   // }
   //  WALBERLA_LOG_INFO( " viscosity before non - D : " << retVal );
   // retVal /= 1e24;
   //  WALBERLA_LOG_INFO( " viscosity after non - D : " << retVal );
   return retVal;
};
class TALASimulation
{
 public:
   TALASimulation( const walberla::Config::BlockHandle& mainConf_,
                   std::shared_ptr< PrimitiveStorage >  storage_,
                   uint_t                               minLevel_,
                   uint_t                               maxLevel_ )
   : mainConf( mainConf_ )
   , storage( storage_ )
   , minLevel( minLevel_ )
   , maxLevel( maxLevel_ )
   , vecMassOperator( storage, minLevel_, maxLevel_ )
   , massOperator( storage, minLevel_, maxLevel_ )
   {
      endTime             = mainConf.getParameter< real_t >( "simulationTime" );
      params.maxTimeSteps = mainConf.getParameter< uint_t >( "maxTimeSteps" );

      params.vtkWriteFrequency = mainConf.getParameter< uint_t >( "vtkWriteFrequency" );

      params.AiniPerturb = mainConf.getParameter< real_t >( "AiniPerturb" );

      params.Ra = mainConf.getParameter< real_t >( "RayleighNumber" );
      params.k_ = mainConf.getParameter< real_t >( "k_" );

      params.minresIter   = mainConf.getParameter< uint_t >( "stokesMinresIter" );
      params.minresRelTol = mainConf.getParameter< real_t >( "stokesMinresTol" );

      params.gmresIter = mainConf.getParameter< uint_t >( "transportGmresIter" );
      params.gmresTol  = mainConf.getParameter< real_t >( "transportGmresTol" );

      tempIni = [this]( const Point3D& x ) {
         return ( 1 - x[2] ) +
                params.AiniPerturb * std::sin( 2.0 * walberla::math::pi * x[0] ) * std::sin( walberla::math::pi * x[2] );
      };

      tempDevBC = []( const Point3D& x ) {
         if ( std::abs( x[2] - 1.0 ) < 1e-6 )
         {
            return 0.0;
         }
         else if ( std::abs( x[2] ) < 1e-6 )
         {
            return 1.0;
         }
         else
         {
            return 0.0;
         }
      };

      rhoFunc = []( const Point3D& x ) { return std::exp( 1.0 - x[2] ); };

      BoundaryCondition bcVelocity, bcPressure, bcVelocityX, bcVelocityY, bcVelocityZ;

      bcTemp.createAllInnerBC();
      bcTemp.createDirichletBC( "DirichletBottomAndTop",
                                { BoundaryMarkers::Top,
                                  BoundaryMarkers::Bottom,
                                  BoundaryMarkers::BLCorner,
                                  BoundaryMarkers::BRCorner,
                                  BoundaryMarkers::BFCorner,
                                  BoundaryMarkers::BBCorner,
                                  BoundaryMarkers::TLCorner,
                                  BoundaryMarkers::TRCorner,
                                  BoundaryMarkers::TFCorner,
                                  BoundaryMarkers::TBCorner,
                                  BoundaryMarkers::PointCorners } );
      bcTemp.createNeumannBC( "NeumannLeftRightFrontBack",
                              { BoundaryMarkers::Left,
                                BoundaryMarkers::Right,
                                BoundaryMarkers::Front,
                                BoundaryMarkers::Back,
                                BoundaryMarkers::FRCorner,
                                BoundaryMarkers::FLCorner,
                                BoundaryMarkers::BaLCorner,
                                BoundaryMarkers::BaRCorner } );

      bcVelocity.createAllInnerBC();
      bcVelocity.createFreeslipBC(
          "AllFreeslip", { BoundaryMarkers::Top, BoundaryMarkers::Bottom, BoundaryMarkers::Left, BoundaryMarkers::Right } );
      bcVelocity.createDirichletBC( "DirichletCorners", { BoundaryMarkers::PointCorners } );

      bcVelocityX.createAllInnerBC();
      bcVelocityX.createDirichletBC( "LRDirichlet",
                                     { BoundaryMarkers::Left,
                                       BoundaryMarkers::Right,
                                       BoundaryMarkers::FRCorner,
                                       BoundaryMarkers::FLCorner,
                                       BoundaryMarkers::BaRCorner,
                                       BoundaryMarkers::BaLCorner,
                                       BoundaryMarkers::TRCorner,
                                       BoundaryMarkers::TLCorner,
                                       BoundaryMarkers::BRCorner,
                                       BoundaryMarkers::BLCorner,
                                       BoundaryMarkers::PointCorners } );
      bcVelocityX.createNeumannBC( "TBFBNeumann",
                                   { BoundaryMarkers::Top,
                                     BoundaryMarkers::Bottom,
                                     BoundaryMarkers::Front,
                                     BoundaryMarkers::Back,
                                     BoundaryMarkers::BFCorner,
                                     BoundaryMarkers::BBCorner,
                                     BoundaryMarkers::TFCorner,
                                     BoundaryMarkers::TBCorner } );

      bcVelocityY.createAllInnerBC();
      bcVelocityY.createDirichletBC( "FBDirichlet",
                                     { BoundaryMarkers::Front,
                                       BoundaryMarkers::Back,
                                       BoundaryMarkers::FRCorner,
                                       BoundaryMarkers::FLCorner,
                                       BoundaryMarkers::BaRCorner,
                                       BoundaryMarkers::BaLCorner,
                                       BoundaryMarkers::BFCorner,
                                       BoundaryMarkers::BBCorner,
                                       BoundaryMarkers::TFCorner,
                                       BoundaryMarkers::TBCorner,
                                       BoundaryMarkers::PointCorners } );
      bcVelocityY.createNeumannBC( "LRTBNeumann",
                                   { BoundaryMarkers::Left,
                                     BoundaryMarkers::Right,
                                     BoundaryMarkers::Top,
                                     BoundaryMarkers::Bottom,
                                     BoundaryMarkers::TRCorner,
                                     BoundaryMarkers::TLCorner,
                                     BoundaryMarkers::BRCorner,
                                     BoundaryMarkers::BLCorner } );

      bcVelocityZ.createAllInnerBC();
      bcVelocityZ.createDirichletBC( "TBDirichlet",
                                     { BoundaryMarkers::Top,
                                       BoundaryMarkers::Bottom,
                                       BoundaryMarkers::TRCorner,
                                       BoundaryMarkers::TLCorner,
                                       BoundaryMarkers::BRCorner,
                                       BoundaryMarkers::BLCorner,
                                       BoundaryMarkers::BFCorner,
                                       BoundaryMarkers::BBCorner,
                                       BoundaryMarkers::TFCorner,
                                       BoundaryMarkers::TBCorner,
                                       BoundaryMarkers::PointCorners } );
      bcVelocityZ.createNeumannBC( "LRFBNeumann",
                                   { BoundaryMarkers::Left,
                                     BoundaryMarkers::Right,
                                     BoundaryMarkers::Front,
                                     BoundaryMarkers::Back,
                                     BoundaryMarkers::FLCorner,
                                     BoundaryMarkers::FRCorner,
                                     BoundaryMarkers::BaLCorner,
                                     BoundaryMarkers::BaRCorner } );

      rhoP2 = std::make_shared< P2Function< real_t > >( "rhoP2", storage_, minLevel_, maxLevel_, bcTemp );

      TP1     = std::make_shared< P1Function< real_t > >( "TP1", storage_, minLevel_, maxLevel_ + 1, bcTemp );
      TRhsP1  = std::make_shared< P1Function< real_t > >( "TRhsP1", storage_, minLevel_, maxLevel_ + 1, bcTemp );
      uP1     = std::make_shared< P1VectorFunction< real_t > >( "uP1", storage_, minLevel_, maxLevel_ + 1, bcVelocity );
      uPrevP1 = std::make_shared< P1VectorFunction< real_t > >( "uPrevP1", storage_, minLevel_, maxLevel_ + 1, bcVelocity );

      T     = std::make_shared< P2Function< real_t > >( "T", storage_, minLevel_, maxLevel_, bcTemp );
      TRhs  = std::make_shared< P2Function< real_t > >( "TRhs", storage_, minLevel_, maxLevel_, bcTemp );
      TPrev = std::make_shared< P2Function< real_t > >( "TPrev", storage_, minLevel_, maxLevel_, bcTemp );

      diffusionTermCoeff = std::make_shared< P2Function< real_t > >( "diffusionTermCoeff", storage_, minLevel_, maxLevel_ );

      u     = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "u", storage_, minLevel_, maxLevel_, bcVelocity );
      uTmp  = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uTmp", storage_, minLevel_, maxLevel_, bcVelocity );
      uPrev = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uPrev", storage_, minLevel_, maxLevel_, bcVelocity );
      uRhsStrong =
          std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhsStrong", storage_, minLevel_, maxLevel_, bcVelocity );
      uRhs = std::make_shared< P2P1TaylorHoodFunction< real_t > >( "uRhs", storage_, minLevel_, maxLevel_, bcVelocity );

      mu    = std::make_shared< P2Function< real_t > >( "viscosity", storage_, minLevel_, maxLevel_, bcTemp );
      muInv = std::make_shared< P1Function< real_t > >( "muInv", storage_, minLevel_, maxLevel_, bcTemp );

      // P2Function< real_t > mu( "viscosity", storage_, minLevel_, maxLevel_, bcTemp ); //mu
      // P1Function< real_t > muInv( "muInv", storage_, minLevel_, maxLevel_, bcTemp );

      u->uvw().setBoundaryCondition( bcVelocityX, 0u );
      u->uvw().setBoundaryCondition( bcVelocityY, 1u );
      u->uvw().setBoundaryCondition( bcVelocityZ, 2u );
      // u->p().setBoundaryCondition( bcPressure );

      uRhs->uvw().setBoundaryCondition( bcVelocityX, 0u );
      uRhs->uvw().setBoundaryCondition( bcVelocityY, 1u );
      uRhs->uvw().setBoundaryCondition( bcVelocityZ, 2u );
      // uRhs->p().setBoundaryCondition( bcPressure );

      rhoP2->interpolate( rhoFunc, maxLevel_, All );

      frozenVelocityOp = std::make_shared< FrozenVelocityOperator >( storage_, minLevel_, maxLevel_, *rhoP2 );

      {
         transportTALAOp = std::make_shared< terraneo::P2TransportOperator >( storage_, minLevel_, maxLevel_ );

         transportTALAOp->setVelocity( u );
         transportTALAOp->setVelocityPrev( uPrev );
         transportTALAOp->setTemperature( T );

         diffusionFunc = [this]( const Point3D& ) { return params.k_; };

         transportTALAOp->setDiffusivityCoeff( std::make_shared< std::function< real_t( const Point3D& ) > >( diffusionFunc ) );

         transportTALAOp->setTALADict( {
             { terraneo::TransportOperatorTermKey::ADVECTION_TERM_WITH_APPLY, false },
             { terraneo::TransportOperatorTermKey::DIFFUSION_TERM, true },
             { terraneo::TransportOperatorTermKey::SHEAR_HEATING_TERM, false },
             { terraneo::TransportOperatorTermKey::ADIABATIC_HEATING_TERM, false },
             { terraneo::TransportOperatorTermKey::INTERNAL_HEATING_TERM, false },
             { terraneo::TransportOperatorTermKey::SUPG_STABILISATION, false },
         } );

         transportTALAOp->initializeOperators();

         transportTALAOp->mmocTransport_->setCautionedEvaluate( mainConf.getParameter< bool >( "cautionedEvaluate" ) );

         transportCGSolver =
             std::make_shared< CGSolver< terraneo::P2TransportOperator > >( storage_, minLevel_, maxLevel_, 30, 1e-9 );
         transportCGSolver->setPrintInfo( params.verbose );
      }

      {
         transportTALAOpP1 = std::make_shared< terraneo::P1TransportOperator >( storage_, minLevel_, maxLevel_ + 1 );

         transportTALAOpP1->setVelocity( u );
         transportTALAOpP1->setVelocityPrev( uPrev );
         transportTALAOpP1->setTemperature( TP1 );

         transportTALAOpP1->setDiffusivityCoeff( std::make_shared< std::function< real_t( const Point3D& ) > >( diffusionFunc ) );

         transportTALAOpP1->setTALADict( {
             { terraneo::TransportOperatorTermKey::ADVECTION_TERM_WITH_APPLY, false },
             { terraneo::TransportOperatorTermKey::DIFFUSION_TERM, true },
             { terraneo::TransportOperatorTermKey::SHEAR_HEATING_TERM, false },
             { terraneo::TransportOperatorTermKey::ADIABATIC_HEATING_TERM, false },
             { terraneo::TransportOperatorTermKey::INTERNAL_HEATING_TERM, false },
             { terraneo::TransportOperatorTermKey::SUPG_STABILISATION, false },
         } );

         transportTALAOpP1->initializeOperators();

         transportCGSolverP1 =
             std::make_shared< CGSolver< terraneo::P1TransportOperator > >( storage_, minLevel_, maxLevel_ + 1, 30, 1e-9 );
         transportCGSolverP1->setPrintInfo( params.verbose );
      }

      real_t allowedRelativeMassDifference = mainConf.getParameter< real_t >( "allowedRelativeMassDifference" );

      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         mu->interpolate( viscFunction, { *T }, level, All );
         muInv->interpolate( viscInv, { *TP1 }, level, All );
      }

      stokesOperator = std::make_shared< StokesOperator_T >( storage_, minLevel_, maxLevel_, *mu );
      stokesOperator->getA().computeInverseDiagonalOperatorValues();

      params.cflMax = mainConf.getParameter< real_t >( "cflMax" );

      schurOperator = std::make_shared< SchurOperator_T >( storage_, minLevel_, maxLevel_, *muInv );

      stokesSolverMG = solvertemplates::fgmresMGSolver< StokesOperator_T, StokesOperator_T::ViscousOperator_T, SchurOperator_T >(
          storage_, minLevel_, maxLevel_, stokesOperator->getA(), *schurOperator );

#if defined( HYTEG_BUILD_WITH_PETSC )
      stokesDirectSolverPetsc = std::make_shared< PETScLUSolver< StokesOperator_T > >( storage_, maxLevel_ );
#endif

      std::string outputFilename = mainConf.getParameter< std::string >( "outputFilename" );
      std::string outputPath     = mainConf.getParameter< std::string >( "outputPath" );

      adios2Output = std::make_shared< AdiosWriter >( outputPath, outputFilename, "", storage );

      adios2Output->add( *u );
      adios2Output->add( *uRhs );
      adios2Output->add( *uTmp );
      adios2Output->add( *T );
      adios2Output->add( *mu );
   }

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > viscFunction =
       [&]( const Point3D& x, const std::vector< real_t >& Temperature ) {
          real_t retVal = viscosityFunction( x, Temperature );
          return retVal;
       };

   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > viscInv =
       [&]( const Point3D& x, const std::vector< real_t >& Temperature ) {
          real_t retVal = 1.0 / viscosityFunction( x, Temperature );
          return retVal;
       };

   void solveU();
   void solveT();
   void step();
   void solve();
   void writeZProfile();
   void writeVTK( uint_t timestep = 0 ) { adios2Output->write( maxLevel, timestep ); }

 private:
   const walberla::Config::BlockHandle& mainConf;

   std::shared_ptr< PrimitiveStorage > storage;
   uint_t                              minLevel, maxLevel;

   std::shared_ptr< P2Function< real_t > > rhoP2;

   BoundaryCondition bcTemp;

   std::shared_ptr< P1Function< real_t > > TP1;
   std::shared_ptr< P1Function< real_t > > TRhsP1;

   std::shared_ptr< P1VectorFunction< real_t > > uP1;
   std::shared_ptr< P1VectorFunction< real_t > > uPrevP1;

   std::shared_ptr< P2Function< real_t > > T;
   std::shared_ptr< P2Function< real_t > > TPrev;
   std::shared_ptr< P2Function< real_t > > TRhs;

   std::shared_ptr< P2Function< real_t > > mu;
   std::shared_ptr< P1Function< real_t > > muInv;

   std::shared_ptr< P2Function< real_t > > diffusionTermCoeff;

   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > u;
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > uTmp;
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > uPrev;
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > uRhsStrong;
   std::shared_ptr< P2P1TaylorHoodFunction< real_t > > uRhs;

   P2ElementwiseBlendingVectorMassOperator vecMassOperator;
   operatorgeneration::P2ElementwiseMass   massOperator;

   std::shared_ptr< StokesOperator_T > stokesOperator;
   std::shared_ptr< SchurOperator_T >  schurOperator;

   std::shared_ptr< terraneo::P2TransportOperator > transportTALAOp;
   std::shared_ptr< terraneo::P1TransportOperator > transportTALAOpP1;

   std::shared_ptr< FrozenVelocityOperator > frozenVelocityOp;

   std::shared_ptr< P1toP1LinearProlongation< real_t > > p1LinearProlongation;

   // Solvers
   std::shared_ptr< Solver< StokesOperator_T > > stokesSolverMG;
   std::shared_ptr< Solver< StokesOperator_T > > stokesDirectSolverPetsc;

   std::shared_ptr< CGSolver< terraneo::P2TransportOperator > > transportCGSolver;
   std::shared_ptr< CGSolver< terraneo::P1TransportOperator > > transportCGSolverP1;

   ParameterContainer params;

   uint_t iTimeStep = 0U;

   real_t simulationTime = 0.0, endTime = 1.0;

   std::function< real_t( const Point3D& ) > diffusionFunc;

   std::function< real_t( const Point3D& ) > tempIni, tempDevBC, rhoFunc, rhoInvFunc, TRefFunc;

   std::function< void( const Point3D&, Point3D& ) > normalsFS;

   // Output

   std::shared_ptr< AdiosWriter > adios2Output;
};

void TALASimulation::solveU()
{
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STARTING STOKES SOLVER" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   uRhsStrong->uvw().component( 0u ).interpolate( 0.0, maxLevel, All );
   uRhsStrong->uvw().component( 1u ).interpolate( 0.0, maxLevel, All );
   uRhsStrong->uvw().component( 2u ).interpolate( params.Ra, maxLevel, All );

   uRhsStrong->uvw().component( 0u ).multElementwise( { uRhsStrong->uvw().component( 0u ), *T }, maxLevel, All );
   uRhsStrong->uvw().component( 1u ).multElementwise( { uRhsStrong->uvw().component( 1u ), *T }, maxLevel, All );
   uRhsStrong->uvw().component( 2u ).multElementwise( { uRhsStrong->uvw().component( 2u ), *T }, maxLevel, All );

   vecMassOperator.apply( uRhsStrong->uvw(), uRhs->uvw(), maxLevel, All );
   // frozenVelocityOp->apply( u->uvw(), uRhs->p(), maxLevel, All );

   for ( uint_t level = minLevel; level <= maxLevel; level++ )
   {
      mu->interpolate( viscFunction, { *T }, level, All );
      muInv->interpolate( viscInv, { *TP1 }, level, All );
   }
   stokesOperator->getA().computeInverseDiagonalOperatorValues();

   if ( mainConf.getParameter< bool >( "stokesDirectSolver" ) )
   {
#if defined( HYTEG_BUILD_WITH_PETSC )
      stokesDirectSolverPetsc->solve( *stokesOperator, *u, *uRhs, maxLevel );
#else
      WALBERLA_ABORT( "PETSc not build but direct solver requested" );
#endif
   }
   else
   {
      stokesSolverMG->solve( *stokesOperator, *u, *uRhs, maxLevel );
   }

   // stokesSolverPetsc->solve( *stokesOperator, *u, *uRhs, maxLevel );
   vertexdof::projectMean( u->p(), maxLevel );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STOKES SOLVER DONE!" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

void TALASimulation::solveT()
{
   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "TRANSPORT SOLVER STARTED!" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   transportTALAOp->calculateTimestep( params.cflMax );
   transportTALAOpP1->setTimestep( transportTALAOp->timestep );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "STARTING TRANSPORT SOLVER with dt = %2.6e", transportTALAOp->timestep ) );

   T->assign( { 1.0 }, { *TPrev }, maxLevel, All );

   if ( mainConf.getParameter< bool >( "p1MMOC" ) )
   {
      P2toP1Conversion( *T, *TP1, maxLevel + 1, All );
      P2toP1Conversion( u->uvw(), *uP1, maxLevel + 1, All );
      P2toP1Conversion( uPrev->uvw(), *uPrevP1, maxLevel + 1, All );

      transportTALAOpP1->mmocTransport_->step(
          *TP1, *uP1, *uPrevP1, maxLevel + 1, hyteg::Inner | hyteg::NeumannBoundary, transportTALAOp->timestep, 1u );

      if ( mainConf.getParameter< bool >( "p1Diffusion" ) )
      {
         TP1->interpolate( tempDevBC, maxLevel + 1, DirichletBoundary );

         transportTALAOpP1->applyRHS( *TRhsP1, maxLevel + 1, All );
         transportCGSolverP1->solve( *transportTALAOpP1, *TP1, *TRhsP1, maxLevel + 1 );

         P1toP2Conversion( *TP1, *T, maxLevel, All );
      }
      else
      {
         P1toP2Conversion( *TP1, *T, maxLevel, All );

         T->interpolate( tempDevBC, maxLevel, DirichletBoundary );

         transportTALAOp->applyRHS( *TRhs, maxLevel, All );
         transportCGSolver->solve( *transportTALAOp, *T, *TRhs, maxLevel );
      }
   }
   else if ( mainConf.getParameter< bool >( "p1Diffusion" ) )
   {
      transportTALAOp->mmocTransport_->step(
          *T, u->uvw(), uPrev->uvw(), maxLevel, hyteg::Inner | hyteg::NeumannBoundary, transportTALAOp->timestep, 1u );

      P2toP1Conversion( *T, *TP1, maxLevel + 1, All );

      TP1->interpolate( tempDevBC, maxLevel + 1, DirichletBoundary );

      transportTALAOpP1->applyRHS( *TRhsP1, maxLevel + 1, All );
      transportCGSolverP1->solve( *transportTALAOpP1, *TP1, *TRhsP1, maxLevel + 1 );

      P1toP2Conversion( *TP1, *T, maxLevel, All );
   }
   else
   {
      transportTALAOp->mmocTransport_->step(
          *T, u->uvw(), uPrev->uvw(), maxLevel, hyteg::Inner | hyteg::NeumannBoundary, transportTALAOp->timestep, 1u );

      T->interpolate( tempDevBC, maxLevel, DirichletBoundary );

      transportTALAOp->applyRHS( *TRhs, maxLevel, All );
      transportCGSolver->solve( *transportTALAOp, *T, *TRhs, maxLevel );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "TRANSPORT SOLVER DONE!" ) );
   WALBERLA_LOG_INFO_ON_ROOT( "" );
}

void TALASimulation::step()
{
   real_t vMax = u->uvw().getMaxComponentMagnitude( maxLevel, All );
   real_t vRMS = std::sqrt( u->uvw().dotGlobal( u->uvw(), maxLevel ) / u->uvw().getNumberOfGlobalDoFs( maxLevel ) );
   real_t hMax = MeshQuality::getMaximalEdgeLength( storage, maxLevel );

   real_t Pe = vLmaxCuboid__ * vRMS / ( params.k_ );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Peclet number (RMS velocity) = %f", Pe ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Maximum velocity = %f", vMax ) );

   uint_t nCouplingIter = mainConf.getParameter< uint_t >( "nCouplingIter" );

   for ( uint_t couplingIter = 0u; couplingIter < nCouplingIter; couplingIter++ )
   {
      solveT();
      solveU();
   }

   transportTALAOp->incrementTimestep();

   TPrev->assign( { 1.0 }, { *T }, maxLevel, All );
   uPrev->assign( { 1.0 }, { *u }, maxLevel, All );

   // WALBERLA_LOG_INFO( "Temperature ?? " );
   // WALBERLA_LOG_INFO( T );
}

void TALASimulation::solve()
{
   T->interpolate( tempIni, maxLevel, Inner | NeumannBoundary );
   T->interpolate( tempDevBC, maxLevel, DirichletBoundary );

   TPrev->assign( { 1.0 }, { *T }, maxLevel, All );

   solveU();

   uPrev->uvw().assign( { 1.0 }, { u->uvw() }, maxLevel, All );

   writeZProfile();
   writeVTK( iTimeStep );

   while ( simulationTime < endTime && iTimeStep < params.maxTimeSteps )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Starting step %d at time = %f!", iTimeStep, simulationTime ) );
      WALBERLA_LOG_INFO_ON_ROOT( "" );

      step();

      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "Step done!" ) );

      iTimeStep++;

      simulationTime += transportTALAOp->timestep;

      if ( iTimeStep % params.vtkWriteFrequency == 0 )
      {
         writeZProfile();
         writeVTK( iTimeStep );
      }
   }
}

void TALASimulation::writeZProfile()
{
   uint_t nz = mainConf.getParameter< uint_t >( "nz" );

   std::string profilePath = mainConf.getParameter< std::string >( "profilePath" );
   std::string fileName    = walberla::format( "%s/Temperature_%d.txt", profilePath.c_str(), iTimeStep );

   auto TProfile =
       terraneo::computeRadialProfile( *T, 0.0, vLmaxCuboid__, nz + 1, maxLevel, []( const Point3D& x ) { return x[2]; } );

   TProfile.logToFile( fileName, "T" );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#if defined( HYTEG_BUILD_WITH_PETSC )
   PETScManager petscManager( &argc, &argv );
#endif

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./CubeConvection.prm" );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   WALBERLA_ROOT_SECTION()
   {
      mainConf.listParameters();
   }

   const uint_t nx = mainConf.getParameter< uint_t >( "nx" );
   const uint_t ny = mainConf.getParameter< uint_t >( "ny" );
   const uint_t nz = mainConf.getParameter< uint_t >( "nz" );

   const uint_t minLevel = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel = mainConf.getParameter< uint_t >( "maxLevel" );

   auto meshInfo = hyteg::MeshInfo::meshCuboid(
       Point3D( 0.0, 0.0, 0.0 ), Point3D( hLmaxCuboid__, hLmaxCuboid__, vLmaxCuboid__ ), nx, ny, nz );

   auto setupStorage = std::make_shared< hyteg::SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::Bottom, bottomMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::Left, leftMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::Right, rightMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::Top, topMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::Front, frontMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::Back, backMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::BRCorner, BRMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::BLCorner, BLMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::BFCorner, BFMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::BBCorner, BBMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::TRCorner, TRMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::TLCorner, TLMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::TFCorner, TFMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::TBCorner, TBMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::FRCorner, FRMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::FLCorner, FLMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::BaRCorner, BaRMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::BaLCorner, BaLMarker );
   setupStorage->setMeshBoundaryFlagsByVertexLocation( BoundaryMarkers::PointCorners, cornersMarker );

   auto storage = std::make_shared< hyteg::PrimitiveStorage >( *setupStorage, 1u );

   uint_t nMacroCells = storage->getNumberOfGlobalCells();

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "Macro Cells = " << nMacroCells );
   WALBERLA_LOG_INFO_ON_ROOT( "" );

   TALASimulation simulation( mainConf, storage, minLevel, maxLevel );

   simulation.solve();

   return 0;
}
