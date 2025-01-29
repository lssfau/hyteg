/*
 * Copyright (c) 2017-2023 Ponsuganth Ilangovan P, Marcus Mohr
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

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/MeshQuality.hpp"
#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/StrongFreeSlipWrapper.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseBlendingFullViscousOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseBlendingStokesOperator.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesProlongation.hpp"
#include "hyteg/gridtransferoperators/P1P1StokesToP1P1StokesRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ProjectNormalOperator.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"
#include "hyteg/python/PythonCallingWrapper.hpp"

#include "mixed_operator/VectorMassOperator.hpp"

using walberla::real_t;

using namespace hyteg;

/***************************/
/** TYPEDEFs on operators **/
/***************************/
typedef P2P1TaylorHoodFunction< real_t > StokesFunction;
typedef P2VectorFunction< real_t >       VelocityFunction;
// typedef P2ElementwiseBlendingFullConstantViscousOperator ViscousVelocityBlockOperator;
// typedef P2P1ElementwiseBlendingStokesOperatorGeneric< ViscousVelocityBlockOperator::ViscousVelocityBlock_0_0,
//                                                       ViscousVelocityBlockOperator >
//                                                          StokesOperator;

typedef P2P1ElementwiseBlendingFullViscousStokesOperator StokesOperator;
// typedef P2P1ElementwiseBlendingStokesOperator StokesOperator;

typedef P2ProjectNormalOperator ProjectionOperator;

typedef P2ElementwiseBlendingVectorMassOperator VelocityVectorMassOperator;
// typedef P2ElementwiseVectorMassOperator   VelocityVectorMassOperator;
typedef P2ElementwiseBlendingMassOperator ScalarMassOperatorVelocity;

typedef P2P1StokesToP2P1StokesProlongation StokesProlongation;

typedef P1Function< real_t >              PressureFunction;
typedef P1ElementwiseBlendingMassOperator ScalarMassOperator;
typedef P2ElementwiseBlendingMassOperator ScalarMassOperatorP2;
/***************************/

// This initializes the Python calling wrapper
PythonCallingWrapper pythonWrapperGlobal( "./",
                                          "analyticalSolutionsSphere",
                                          { "getDirichletVelocitySmooth2d",
                                            "getDirichletPressureSmooth2d",
                                            "getFreeslipVelocitySmooth2d",
                                            "getFreeslipPressureSmooth2d",
                                            "getFreeZeroslipVelocitySmooth2d",
                                            "getFreeZeroslipPressureSmooth2d",
                                            "getDirichletVelocityDelta2d",
                                            "getDirichletPressureDelta2d",
                                            "getFreeslipVelocityDelta2d",
                                            "getFreeslipPressureDelta2d",
                                            "getFreeZeroslipVelocityDelta2d",
                                            "getFreeZeroslipPressureDelta2d",
                                            "getDeltaRho2d" } );

template < typename FunctionType, typename MassOperator >
real_t normL2( const FunctionType& u, const FunctionType& tmp, const MassOperator& M, uint_t level, DoFType flag )
{
   M.apply( u, tmp, level, flag );
   return std::sqrt( u.dotGlobal( tmp, level, flag ) );
}

enum BoundaryConditionType
{
   ALL_DIRICHLET,
   ALL_FREESLIP,
   MIXED_DIRICHLET_AND_FREESLIP
};

template < bool deltaForcing, BoundaryConditionType boundaryConditionType, int component >
real_t uvwSolution( const hyteg::Point3D& x )
{
   if ( deltaForcing )
   {
      if ( boundaryConditionType == ALL_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeslipVelocityDelta2d" )[component];
      }
      else if ( boundaryConditionType == MIXED_DIRICHLET_AND_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeZeroslipVelocityDelta2d" )[component];
      }
      else
      {
         return pythonWrapperGlobal.getParameter( x, "getDirichletVelocityDelta2d" )[component];
      }
   }
   else
   {
      if ( boundaryConditionType == ALL_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeslipVelocitySmooth2d" )[component];
      }
      else if ( boundaryConditionType == MIXED_DIRICHLET_AND_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeZeroslipVelocitySmooth2d" )[component];
      }
      else
      {
         return pythonWrapperGlobal.getParameter( x, "getDirichletVelocitySmooth2d" )[component];
      }
   }
}

template < bool deltaForcing, BoundaryConditionType boundaryConditionType >
real_t pSolution( const hyteg::Point3D& x )
{
   if ( deltaForcing )
   {
      if ( boundaryConditionType == ALL_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeslipPressureDelta2d" )[0];
      }
      else if ( boundaryConditionType == MIXED_DIRICHLET_AND_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeZeroslipPressureDelta2d" )[0];
      }
      else
      {
         return pythonWrapperGlobal.getParameter( x, "getDirichletPressureDelta2d" )[0];
      }
   }
   else
   {
      if ( boundaryConditionType == ALL_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeslipPressureSmooth2d" )[0];
      }
      else if ( boundaryConditionType == MIXED_DIRICHLET_AND_FREESLIP )
      {
         return pythonWrapperGlobal.getParameter( x, "getFreeZeroslipPressureSmooth2d" )[0];
      }
      else
      {
         return pythonWrapperGlobal.getParameter( x, "getDirichletPressureSmooth2d" )[0];
      }
   }
}

real_t diracDelta( const hyteg::Point3D& x, real_t rDash )
{
   if ( std::abs( std::sqrt( x[0] * x[0] + x[1] * x[1] ) - rDash ) < 1e-4 )
      return 1.0;
   else
      return 0.0;
}

template < bool deltaForcing, int component >
real_t forcingFunction( const hyteg::Point3D& x )
{
   real_t r     = std::sqrt( x[0] * x[0] + x[1] * x[1] );
   real_t theta = std::atan2( x[1], x[0] );
   // real_t R_plus = 2.22;
   real_t rDash = 1.72;

   if ( deltaForcing )
   {
      real_t rhoDash = std::cos( 2 * theta );
      return -diracDelta( x, rDash ) * rhoDash * x[component] / r;
   }
   else
   {
      real_t rhoDash = pythonWrapperGlobal.getParameter( x, "getDeltaRho2d" )[0];
      return -rhoDash * x[component] / r;
   }
}

void removeRotationalModes( P2ElementwiseBlendingMassOperator& massOperator,
                            P2VectorFunction< real_t >&        u,
                            P2VectorFunction< real_t >&        rtheta,
                            P2Function< real_t >&              temp,
                            uint_t                             level )
{
   if ( temp.getStorage()->hasGlobalCells() )
   {
      std::function< real_t( const Point3D& ) > rSquareFunc = []( const Point3D& x ) { return x[0] * x[0] + x[1] * x[1]; };

      rtheta[2].interpolate( rSquareFunc, level, All );
      massOperator.apply( rtheta[2], temp, level, All );
      real_t rSquare = temp.dotGlobal( rtheta[2], level, All );

      /***************************************************************************/
      // Z axis mode

      std::function< real_t( const Point3D& ) > zAxisModeX = []( const Point3D& x ) { return -x[1]; };
      std::function< real_t( const Point3D& ) > zAxisModeY = []( const Point3D& x ) { return x[0]; };

      rtheta[0].interpolate( zAxisModeX, level, All );
      rtheta[1].interpolate( zAxisModeY, level, All );
      rtheta[2].interpolate( 0.0, level, All );

      rtheta[0].multElementwise( { rtheta[0], u[0] }, level, All );
      rtheta[1].multElementwise( { rtheta[1], u[1] }, level, All );

      rtheta[2].assign( { 1.0, 1.0 }, { rtheta[0], rtheta[1] }, level, All );
      massOperator.apply( rtheta[2], temp, level, All );
      real_t rThetaDotUZ = temp.dotGlobal( rtheta[2], level, All );

      rtheta[0].interpolate( zAxisModeX, level, All );
      rtheta[1].interpolate( zAxisModeY, level, All );
      rtheta[2].interpolate( 0.0, level, All );

      u.assign( { 1.0, -rThetaDotUZ / rSquare }, { u, rtheta }, level, All );

      /***************************************************************************/
      // X axis mode

      std::function< real_t( const Point3D& ) > xAxisModeX = []( const Point3D& x ) { return -x[2]; };
      std::function< real_t( const Point3D& ) > xAxisModeY = []( const Point3D& x ) { return x[1]; };

      rtheta[1].interpolate( xAxisModeX, level, All );
      rtheta[2].interpolate( xAxisModeY, level, All );
      rtheta[0].interpolate( 0.0, level, All );

      rtheta[1].multElementwise( { rtheta[1], u[1] }, level, All );
      rtheta[2].multElementwise( { rtheta[2], u[2] }, level, All );

      rtheta[0].assign( { 1.0, 1.0 }, { rtheta[1], rtheta[2] }, level, All );
      massOperator.apply( rtheta[0], temp, level, All );
      real_t rThetaDotUX = temp.dotGlobal( rtheta[0], level, All );

      rtheta[1].interpolate( xAxisModeX, level, All );
      rtheta[2].interpolate( xAxisModeY, level, All );
      rtheta[0].interpolate( 0.0, level, All );

      u.assign( { 1.0, -rThetaDotUX / rSquare }, { u, rtheta }, level, All );

      /***************************************************************************/
      // Y axis mode

      std::function< real_t( const Point3D& ) > yAxisModeX = []( const Point3D& x ) { return -x[0]; };
      std::function< real_t( const Point3D& ) > yAxisModeY = []( const Point3D& x ) { return x[2]; };

      rtheta[2].interpolate( yAxisModeX, level, All );
      rtheta[0].interpolate( yAxisModeY, level, All );
      rtheta[1].interpolate( 0.0, level, All );

      rtheta[2].multElementwise( { rtheta[2], u[2] }, level, All );
      rtheta[0].multElementwise( { rtheta[0], u[0] }, level, All );

      rtheta[1].assign( { 1.0, 1.0 }, { rtheta[2], rtheta[0] }, level, All );
      massOperator.apply( rtheta[1], temp, level, All );
      real_t rThetaDotUY = temp.dotGlobal( rtheta[1], level, All );

      rtheta[2].interpolate( yAxisModeX, level, All );
      rtheta[0].interpolate( yAxisModeY, level, All );
      rtheta[1].interpolate( 0.0, level, All );

      u.assign( { 1.0, -rThetaDotUY / rSquare }, { u, rtheta }, level, All );
   }
   else
   {
      std::function< real_t( const Point3D& ) > rSquareFunc = []( const Point3D& x ) { return x[0] * x[0] + x[1] * x[1]; };

      rtheta[1].interpolate( rSquareFunc, level, All );
      massOperator.apply( rtheta[1], temp, level, All );
      real_t rSquare = temp.dotGlobal( rtheta[1], level, All );

      /***************************************************************************/
      // Z axis mode

      std::function< real_t( const Point3D& ) > zAxisModeX = []( const Point3D& x ) { return -x[1]; };
      std::function< real_t( const Point3D& ) > zAxisModeY = []( const Point3D& x ) { return x[0]; };

      rtheta[0].interpolate( zAxisModeX, level, All );
      rtheta[1].interpolate( zAxisModeY, level, All );

      rtheta[0].multElementwise( { rtheta[0], u[0] }, level, All );
      rtheta[1].multElementwise( { rtheta[1], u[1] }, level, All );

      rtheta[1].assign( { 1.0, 1.0 }, { rtheta[0], rtheta[1] }, level, All );
      massOperator.apply( rtheta[1], temp, level, All );
      real_t rThetaDotUZ = temp.dotGlobal( rtheta[1], level, All );

      rtheta[0].interpolate( zAxisModeX, level, All );
      rtheta[1].interpolate( zAxisModeY, level, All );

      u.assign( { 1.0, -rThetaDotUZ / rSquare }, { u, rtheta }, level, All );
   }
}

std::array< real_t, 3 >
    convAnalysis( const walberla::Config::BlockHandle& mainConf, uint_t minLevel, uint_t maxLevel, real_t& hMax )
{
   // const real_t rMin = walberla::math::pi, rMax = 2 * walberla::math::pi;
   const real_t rMin = 1.22, rMax = 2.22;

   real_t rMean = ( rMin + rMax ) / real_c( 2.0 );

   // real_t rDp = 1.97;
   // real_t rDm = 1.47;

   uint_t nTan = mainConf.getParameter< uint_t >( "annulusNTan" );
   uint_t nRad = mainConf.getParameter< uint_t >( "annulusNRad" );

   MeshInfo meshInfo = MeshInfo::meshAnnulus( rMin, rMax, MeshInfo::CROSS, nTan, nRad );

   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   AnnulusMap::setMap( *setupStorage );

   std::function< bool( const Point3D& ) > faceMarker = [rMean]( const Point3D& x ) {
      if ( std::abs( std::sqrt( x[0] * x[0] + x[1] * x[1] ) - rMean ) < 1e-6 )
      {
         return true;
      }
      else
      {
         return false;
      }
      // if ( std::abs( std::sqrt( x[0] * x[0] + x[1] * x[1] ) - rDp ) < 1e-6 ||
      //      std::abs( std::sqrt( x[0] * x[0] + x[1] * x[1] ) - rDm ) < 1e-6 )
      // {
      //    return true;
      // }
      // else
      // {
      //    return false;
      // }
   };

   std::function< bool( const Point3D& ) > faceMarkerInner = [rMin]( const Point3D& x ) {
      if ( std::abs( std::sqrt( x[0] * x[0] + x[1] * x[1] ) - rMin ) < 1e-6 )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   std::function< bool( const Point3D& ) > faceMarkerOuter = [rMax]( const Point3D& x ) {
      if ( std::abs( std::sqrt( x[0] * x[0] + x[1] * x[1] ) - rMax ) < 1e-6 )
      {
         return true;
      }
      else
      {
         return false;
      }
   };

   enum boundaryConditionFlagsLocal
   {
      INNER_SPECIAL = 17
   };

   setupStorage->setMeshBoundaryFlagsByCentroidLocation( INNER_SPECIAL, faceMarker );

   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage, 3 );

   BoundaryCondition bcVelocity;

   bcVelocity.createAllInnerBC();

   if ( mainConf.getParameter< bool >( "freeslip" ) )
   {
      bcVelocity.createFreeslipBC( "FreeslipInner", { MeshInfo::hollowFlag::flagInnerBoundary } );
      bcVelocity.createFreeslipBC( "FreeslipOuter", { MeshInfo::hollowFlag::flagOuterBoundary } );
   }
   else if ( mainConf.getParameter< bool >( "mixed" ) )
   {
      bcVelocity.createDirichletBC( "DirichletOuter", { MeshInfo::hollowFlag::flagOuterBoundary } );
      bcVelocity.createFreeslipBC( "FreeslipInner", { MeshInfo::hollowFlag::flagInnerBoundary } );
   }
   else
   {
      bcVelocity.createDirichletBC( "DirichletInner", { MeshInfo::hollowFlag::flagInnerBoundary } );
      bcVelocity.createDirichletBC( "DirichletOuter", { MeshInfo::hollowFlag::flagOuterBoundary } );
   }

   hMax = MeshQuality::getMaximalEdgeLength( storage, maxLevel );
   // real_t hMax = MeshQuality::getMaximalEdgeLength( storage, maxLevel );

   // real_t hMean = ( hMin + hMax ) / 2;

   // VTKOutput vtkOutput( "./output", "2DAnnulus", storage );

   std::shared_ptr< VTKOutput > vtkOutput;
   std::shared_ptr< VTKOutput > vtkOutputHigher;

   if ( mainConf.getParameter< bool >( "delta" ) )
   {
      if ( mainConf.getParameter< bool >( "freeslip" ) )
      {
         vtkOutput       = std::make_shared< hyteg::VTKOutput >( "./output", "2DAnnulusDeltaFS", storage );
         vtkOutputHigher = std::make_shared< hyteg::VTKOutput >( "./output", "2DAnnulusDeltaFSError", storage );
      }
      else if ( mainConf.getParameter< bool >( "mixed" ) )
      {
         vtkOutput       = std::make_shared< hyteg::VTKOutput >( "./output", "2DAnnulusDeltaMX", storage );
         vtkOutputHigher = std::make_shared< hyteg::VTKOutput >( "./output", "2DAnnulusDeltaMXError", storage );
      }
      else
      {
         vtkOutput       = std::make_shared< hyteg::VTKOutput >( "./output", "2DAnnulusDeltaZS", storage );
         vtkOutputHigher = std::make_shared< hyteg::VTKOutput >( "./output", "2DAnnulusDeltaZSError", storage );
      }
   }
   else
   {
      if ( mainConf.getParameter< bool >( "freeslip" ) )
      {
         vtkOutput       = std::make_shared< hyteg::VTKOutput >( "./output", "2DAnnulusSmoothFS", storage );
         vtkOutputHigher = std::make_shared< hyteg::VTKOutput >( "./output", "2DAnnulusSmoothFSError", storage );
      }
      else if ( mainConf.getParameter< bool >( "mixed" ) )
      {
         vtkOutput       = std::make_shared< hyteg::VTKOutput >( "./output", "2DAnnulusSmoothMX", storage );
         vtkOutputHigher = std::make_shared< hyteg::VTKOutput >( "./output", "2DAnnulusSmoothMXError", storage );
      }
      else
      {
         vtkOutput       = std::make_shared< hyteg::VTKOutput >( "./output", "2DAnnulusSmoothZS", storage );
         vtkOutputHigher = std::make_shared< hyteg::VTKOutput >( "./output", "2DAnnulusSmoothZSError", storage );
      }
   }

   StokesFunction   u( "u", storage, minLevel, maxLevel + 1, bcVelocity );
   StokesFunction   rhs( "rhs", storage, minLevel, maxLevel, bcVelocity );
   StokesFunction   rhsStrong( "rhsStrong", storage, minLevel, maxLevel, bcVelocity );
   StokesFunction   residual( "residual", storage, maxLevel, maxLevel, bcVelocity );
   StokesFunction   Au( "Au", storage, maxLevel, maxLevel );
   StokesFunction   AuAnalytical( "AuAnalytical", storage, maxLevel + 1, maxLevel + 1 );
   StokesFunction   uAnalytical( "uAnalytical", storage, maxLevel + 1, maxLevel + 1, bcVelocity );
   StokesFunction   error( "error", storage, minLevel, maxLevel + 1, bcVelocity );
   VelocityFunction tmp( "tmp", storage, minLevel, maxLevel + 1, bcVelocity );
   VelocityFunction Id( "Id", storage, minLevel, maxLevel + 1, bcVelocity );
   PressureFunction tmpP( "tmpP", storage, minLevel, maxLevel + 1 );

   VelocityFunction     rtheta( "u", storage, minLevel, maxLevel, bcVelocity );
   P2Function< real_t > temp( "temp", storage, minLevel, maxLevel );

   vtkOutput->add( u );
   vtkOutput->add( rhs );
   // vtkOutput->add( rhsSurf );

   // real_t uInner = real_c( 4.0 );
   // real_t uOuter = real_c( -8.0 );

   // u.uvw().interpolate( 0, maxLevel, DirichletBoundary );
   // u.p().interpolate( analyticalP, maxLevel, All );

   std::shared_ptr< StokesOperator > stokesOperator = std::make_shared< StokesOperator >( storage, minLevel, maxLevel + 1 );

   auto normalsAnnulus = [rMean]( const Point3D& x, Point3D& n ) {
      real_t r = std::sqrt( x[0] * x[0] + x[1] * x[1] );

      if ( r < rMean )
      {
         n[0] = -x[0] / r;
         n[1] = -x[1] / r;
      }
      else if ( r > rMean )
      {
         n[0] = x[0] / r;
         n[1] = x[1] / r;
      }
   };

   using StokesOperatorFS = hyteg::StrongFreeSlipWrapper< StokesOperator, ProjectionOperator, true >;

   std::shared_ptr< ProjectionOperator > projectNormal =
       std::make_shared< ProjectionOperator >( storage, minLevel, maxLevel, normalsAnnulus );

   real_t minresTol  = mainConf.getParameter< real_t >( "tol" );
   uint_t minresIter = mainConf.getParameter< uint_t >( "minreIter" );

   /*****************************************/
   /****** TRY A DIFFERENT SOLVER HERE ******/
   /*****************************************/
   auto stokesSolver =
       hyteg::solvertemplates::stokesMinResSolver< StokesOperatorFS >( storage, maxLevel, minresTol, minresIter, true );
   /*****************************************/

   StokesOperatorFS stokesOperatorFS( stokesOperator, projectNormal, FreeslipBoundary );

   P2ElementwiseMassOperator  MLineIntegral( storage, minLevel, maxLevel );
   VelocityVectorMassOperator M( storage, minLevel, maxLevel + 1 );
   ScalarMassOperator         scalarM( storage, minLevel, maxLevel + 1 );
   ScalarMassOperatorP2       scalarMP2( storage, minLevel, maxLevel + 1 );

   if ( mainConf.getParameter< bool >( "delta" ) )
   {
      if ( mainConf.getParameter< bool >( "freeslip" ) )
      {
         std::vector< std::function< real_t( const Point3D& ) > > uvwSolutionVec = { uvwSolution< true, ALL_FREESLIP, 0 >,
                                                                                     uvwSolution< true, ALL_FREESLIP, 1 > };
         u.uvw().interpolate( uvwSolutionVec, maxLevel, DirichletBoundary );

         uAnalytical.uvw().interpolate( uvwSolutionVec, maxLevel + 1, All );
         uAnalytical.p().interpolate( pSolution< true, ALL_FREESLIP >, maxLevel + 1, All );
      }
      else if ( mainConf.getParameter< bool >( "mixed" ) )
      {
         std::vector< std::function< real_t( const Point3D& ) > > uvwSolutionVec = {
             uvwSolution< true, MIXED_DIRICHLET_AND_FREESLIP, 0 >, uvwSolution< true, MIXED_DIRICHLET_AND_FREESLIP, 1 > };

         u.uvw().interpolate( uvwSolutionVec, maxLevel, DirichletBoundary );

         uAnalytical.uvw().interpolate( uvwSolutionVec, maxLevel + 1, All );
         uAnalytical.p().interpolate( pSolution< true, MIXED_DIRICHLET_AND_FREESLIP >, maxLevel + 1, All );
      }
      else
      {
         std::vector< std::function< real_t( const Point3D& ) > > uvwSolutionVec = { uvwSolution< true, ALL_DIRICHLET, 0 >,
                                                                                     uvwSolution< true, ALL_DIRICHLET, 1 > };

         u.uvw().interpolate( uvwSolutionVec, maxLevel, DirichletBoundary );

         uAnalytical.uvw().interpolate( uvwSolutionVec, maxLevel + 1, All );
         uAnalytical.p().interpolate( pSolution< true, ALL_DIRICHLET >, maxLevel + 1, All );
      }

      std::vector< std::function< real_t( const Point3D& ) > > rhsVecFunc = { forcingFunction< true, 0 >,
                                                                              forcingFunction< true, 1 > };

      rhsStrong.uvw().interpolate( rhsVecFunc, maxLevel, All );
      real_t scalingValue = real_c( nRad ) * std::pow( 2.0, maxLevel ) * 3.0;
      rhsStrong.uvw().assign( { scalingValue }, { rhsStrong.uvw() }, maxLevel, All );
   }
   else
   {
      if ( mainConf.getParameter< bool >( "freeslip" ) )
      {
         std::vector< std::function< real_t( const Point3D& ) > > uvwSolutionVec = { uvwSolution< false, ALL_FREESLIP, 0 >,
                                                                                     uvwSolution< false, ALL_FREESLIP, 1 > };

         u.uvw().interpolate( uvwSolutionVec, maxLevel, DirichletBoundary );

         uAnalytical.uvw().interpolate( uvwSolutionVec, maxLevel + 1, All );
         uAnalytical.p().interpolate( pSolution< false, ALL_FREESLIP >, maxLevel + 1, All );
      }
      else if ( mainConf.getParameter< bool >( "mixed" ) )
      {
         std::vector< std::function< real_t( const Point3D& ) > > uvwSolutionVec = {
             uvwSolution< false, MIXED_DIRICHLET_AND_FREESLIP, 0 >, uvwSolution< false, MIXED_DIRICHLET_AND_FREESLIP, 1 > };

         u.uvw().interpolate( uvwSolutionVec, maxLevel, DirichletBoundary );

         uAnalytical.uvw().interpolate( uvwSolutionVec, maxLevel + 1, All );
         uAnalytical.p().interpolate( pSolution< false, MIXED_DIRICHLET_AND_FREESLIP >, maxLevel + 1, All );
      }
      else
      {
         std::vector< std::function< real_t( const Point3D& ) > > uvwSolutionVec = { uvwSolution< false, ALL_DIRICHLET, 0 >,
                                                                                     uvwSolution< false, ALL_DIRICHLET, 1 > };

         u.uvw().interpolate( uvwSolutionVec, maxLevel, DirichletBoundary );

         uAnalytical.uvw().interpolate( uvwSolutionVec, maxLevel + 1, All );
         uAnalytical.p().interpolate( pSolution< false, ALL_DIRICHLET >, maxLevel + 1, All );
      }

      std::vector< std::function< real_t( const Point3D& ) > > rhsVecFunc = { forcingFunction< false, 0 >,
                                                                              forcingFunction< false, 1 > };

      rhsStrong.uvw().interpolate( rhsVecFunc, maxLevel, All );
   }

   std::function< real_t( const Point3D& ) > f2DX = []( const Point3D& x ) {
      real_t r     = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      real_t theta = std::atan2( x[1], x[0] );
      return -std::cos( 2 * theta ) * x[0] / ( 2.0 * r );
   };

   std::function< real_t( const Point3D& ) > f2DY = []( const Point3D& x ) {
      real_t r     = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      real_t theta = std::atan2( x[1], x[0] );
      return -std::cos( 2 * theta ) * x[1] / ( 2.0 * r );
   };

   if ( mainConf.getParameter< bool >( "delta" ) )
   {
      if ( mainConf.getParameter< bool >( "surfIntBlend" ) )
      {
         // P2ElementwiseBlendingSurfaceDeltaOperator testOperatorX( storage, maxLevel, maxLevel, f2DX );
         // P2ElementwiseBlendingSurfaceDeltaOperator testOperatorY( storage, maxLevel, maxLevel, f2DY );

         // testOperatorX.applySurface( rhsSurf[0], maxLevel, InnerSpecial );
         // testOperatorY.applySurface( rhsSurf[1], maxLevel, InnerSpecial );

         // rhs.uvw().assign( { 1.0 }, { rhsSurf }, maxLevel, All );

         WALBERLA_ABORT( "Cannot use surface integrals, atleast not in this branch" );
      }
      else if ( mainConf.getParameter< bool >( "surfInt" ) )
      {
         // P2ElementwiseSurfaceDeltaOperator testOperatorX( storage, maxLevel, maxLevel, f2DX );
         // P2ElementwiseSurfaceDeltaOperator testOperatorY( storage, maxLevel, maxLevel, f2DY );

         // testOperatorX.applySurface( rhsSurf[0], maxLevel, InnerSpecial );
         // testOperatorY.applySurface( rhsSurf[1], maxLevel, InnerSpecial );

         // rhs.uvw().assign( { 1.0 }, { rhsSurf }, maxLevel, All );

         WALBERLA_ABORT( "Cannot use surface integrals, atleast not in this branch" );
      }
      else
      {
         M.apply( rhsStrong.uvw(), rhs.uvw(), maxLevel, All );
      }

      // rhs.uvw().assign( { 0.5 }, { rhs.uvw() }, maxLevel, All );
      // M.apply( rhsStrong.uvw(), rhs.uvw(), maxLevel, All );
   }
   else
   {
      M.apply( rhsStrong.uvw(), rhs.uvw(), maxLevel, All );
   }

   // auto directSolver = PETScLUSolver< StokesOperatorFS >( storage, maxLevel );

   uint_t solverType = mainConf.getParameter< uint_t >( "solverType" );
   if ( solverType == 0U )
   {
      // projectNormal->project( rhs, maxLevel, FreeslipBoundary, RotationNormal );
      // directSolver.solve( stokesOperatorFS, u, rhs, maxLevel );
      // projectNormal->project( u, maxLevel, FreeslipBoundary, RotationNormal );

      // vertexdof::projectMean( u.p(), maxLevel );

      WALBERLA_ABORT( "Direct solver will not work, atleast in this branch, as rotation matrix stuff is not merged" );
   }
   else if ( solverType == 1U )
   {
      projectNormal->project( rhs, maxLevel, FreeslipBoundary );
      stokesSolver->solve( stokesOperatorFS, u, rhs, maxLevel );

      vertexdof::projectMean( u.p(), maxLevel );
   }
   // else if ( solverType == 3U ) {}

   /* REMOVE ROTATIONAL MODES */

   // removeRotationalModes( scalarMP2, u.uvw(), rtheta, temp, maxLevel );

   /* ERROR CALCULATION */

   stokesOperator->apply( u, Au, maxLevel, All );

   residual.assign( { real_c( 1.0 ), real_c( -1.0 ) }, { rhs, Au }, maxLevel, All );

   vtkOutput->add( residual );

   stokesOperator->apply( uAnalytical, AuAnalytical, maxLevel + 1, Inner );

   vtkOutput->add( Au );
   vtkOutputHigher->add( AuAnalytical );
   vtkOutput->add( rhsStrong );

   StokesProlongation stokesProlongation;

   stokesProlongation.prolongate( u, maxLevel, All );

   error.assign( { real_c( 1.0 ), real_c( -1.0 ) }, { u, uAnalytical }, maxLevel + 1 );

   vtkOutputHigher->add( uAnalytical );
   vtkOutput->add( error );

   vtkOutputHigher->add( error );

   if ( mainConf.getParameter< bool >( "writeVTK" ) )
   {
      vtkOutput->write( maxLevel );
      vtkOutputHigher->write( maxLevel + 1 );
   }

   real_t errorUV   = normL2( error.uvw(), tmp, M, maxLevel + 1, All );
   real_t errorSurf = normL2( error.uvw(), tmp, M, maxLevel + 1, FreeslipBoundary );
   real_t errorP    = normL2( error.p(), tmpP, scalarM, maxLevel + 1, All );

   // real_t uNorm     = normL2( u.uvw(), rhsStrong.uvw(), M, maxLevel, All );
   // real_t uSurfNorm = normL2( u.uvw(), rhsStrong.uvw(), M, maxLevel, FreeslipBoundary );
   // real_t pNorm     = normL2( u.p(), rhsStrong.p(), scalarM, maxLevel, All );

   real_t residualNorm = std::sqrt( residual.uvw().dotGlobal( residual.uvw(), maxLevel, All ) );

   WALBERLA_ROOT_SECTION() { std::cout << "Residual = " << residualNorm << std::endl << "errorUV = " << errorUV << std::endl; }

   return { errorUV, errorSurf, errorP };
   // return { errorUV / uNorm, errorSurf / uSurfNorm, errorP / pNorm };
}

int main( int argc, char* argv[] )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      cfg->readParameterFile( "./2DAnnulus.prm" );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   WALBERLA_ROOT_SECTION() { mainConf.listParameters(); }

   real_t hMax = real_c( 0.0 );

   const uint_t level = mainConf.getParameter< uint_t >( "level" );

   std::array< real_t, 3 > l2Error = convAnalysis( mainConf, level, level, hMax );

   std::string outputPath = mainConf.getParameter< std::string >( "outputPath" );

   WALBERLA_ROOT_SECTION()
   {
      std::cout << walberla::format( "hMax = %10.10e, errorU = %10.10e", hMax, l2Error[0] ) << std::endl;
      std::cout << walberla::format( "hMax = %10.10e, errorU_Surface = %10.10e", hMax, l2Error[1] ) << std::endl;
      std::cout << walberla::format( "hMax = %10.10e, errorP = %10.10e", hMax, l2Error[2] ) << std::endl;

      if ( mainConf.getParameter< bool >( "calcAndWriteError" ) )
      {
         std::ofstream fileU;
         std::ofstream fileUSurface;
         std::ofstream fileP;

         if ( mainConf.getParameter< bool >( "delta" ) )
         {
            if ( mainConf.getParameter< bool >( "freeslip" ) )
            {
               fileU        = std::ofstream( outputPath + "error_analysis_2D_U_Delta_FS.txt", std::ofstream::app );
               fileUSurface = std::ofstream( outputPath + "error_analysis_2D_USurface_Delta_FS.txt", std::ofstream::app );
               fileP        = std::ofstream( outputPath + "error_analysis_2D_P_Delta_FS.txt", std::ofstream::app );
            }
            else if ( mainConf.getParameter< bool >( "mixed" ) )
            {
               fileU        = std::ofstream( outputPath + "error_analysis_2D_U_Delta_MX.txt", std::ofstream::app );
               fileUSurface = std::ofstream( outputPath + "error_analysis_2D_USurface_Delta_MX.txt", std::ofstream::app );
               fileP        = std::ofstream( outputPath + "error_analysis_2D_P_Delta_MX.txt", std::ofstream::app );
            }
            else
            {
               fileU        = std::ofstream( outputPath + "error_analysis_2D_U_Delta_ZS.txt", std::ofstream::app );
               fileUSurface = std::ofstream( outputPath + "error_analysis_2D_USurface_Delta_ZS.txt", std::ofstream::app );
               fileP        = std::ofstream( outputPath + "error_analysis_2D_P_Delta_ZS.txt", std::ofstream::app );
            }
         }
         else
         {
            if ( mainConf.getParameter< bool >( "freeslip" ) )
            {
               fileU        = std::ofstream( outputPath + "error_analysis_2D_U_Smooth_FS.txt", std::ofstream::app );
               fileUSurface = std::ofstream( outputPath + "error_analysis_2D_USurface_Smooth_FS.txt", std::ofstream::app );
               fileP        = std::ofstream( outputPath + "error_analysis_2D_P_Smooth_FS.txt", std::ofstream::app );
            }
            else if ( mainConf.getParameter< bool >( "mixed" ) )
            {
               fileU        = std::ofstream( outputPath + "error_analysis_2D_U_Smooth_MX.txt", std::ofstream::app );
               fileUSurface = std::ofstream( outputPath + "error_analysis_2D_USurface_Smooth_MX.txt", std::ofstream::app );
               fileP        = std::ofstream( outputPath + "error_analysis_2D_P_Smooth_MX.txt", std::ofstream::app );
            }
            else
            {
               fileU        = std::ofstream( outputPath + "error_analysis_2D_U_Smooth_ZS.txt", std::ofstream::app );
               fileUSurface = std::ofstream( outputPath + "error_analysis_2D_USurface_Smooth_ZS.txt", std::ofstream::app );
               fileP        = std::ofstream( outputPath + "error_analysis_2D_P_Smooth_ZS.txt", std::ofstream::app );
            }
         }

         WALBERLA_ROOT_SECTION()
         {
            fileU << walberla::format( "%10.10e, %10.10e", hMax, l2Error[0] ) << std::endl;
            fileUSurface << walberla::format( "%10.10e, %10.10e", hMax, l2Error[1] ) << std::endl;
            fileP << walberla::format( "%10.10e, %10.10e", hMax, l2Error[2] ) << std::endl;
         }

         fileU.close();
         fileUSurface.close();
         fileP.close();
      }
   }

   return 0;
}
