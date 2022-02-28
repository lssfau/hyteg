/*
 * Copyright (c) 2017-2021 Nils Kohl.
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

#include <cmath>

#include "core/DataTypes.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1EpsilonStokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1P1ElementwiseAffineEpsilonStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseAffineEpsilonStokesOperator.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseConstantCoefficientStokesOperator.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScExportOperatorMatrix.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/Solver.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesBlockDiagonalPreconditioner.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/solvertemplates/StokesSolverTemplates.hpp"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

std::function< real_t( const hyteg::Point3D& ) > plumeExactU = []( const hyteg::Point3D& x ) {
   return std::sin( 2 * pi * x[0] ) * std::cos( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > plumeBCU = []( const hyteg::Point3D& x ) {
   return std::sin( 2 * pi * x[0] ) * std::cos( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > plumeExactV = []( const hyteg::Point3D& x ) {
   return -2.0 * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > plumeBCV = []( const hyteg::Point3D& x ) {
   return -2.0 * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > plumeExactP = []( const hyteg::Point3D& x ) {
   return 2.5 * pi * std::cos( 2 * pi * x[0] ) * std::cos( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > plumeRhsU = []( const hyteg::Point3D& ) { return 0; };
std::function< real_t( const hyteg::Point3D& ) > plumeRhsV = []( const hyteg::Point3D& x ) {
   return -12.5 * pi * pi * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
};

std::function< real_t( const hyteg::Point3D& ) > shellExactU = []( const hyteg::Point3D& x ) {
   return -4 * std::cos( 4 * x[2] );
};
std::function< real_t( const hyteg::Point3D& ) > shellExactV = []( const hyteg::Point3D& x ) { return 8 * std::cos( 8 * x[0] ); };
std::function< real_t( const hyteg::Point3D& ) > shellExactW = []( const hyteg::Point3D& x ) {
   return -2 * std::cos( 2 * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > shellExactP = []( const hyteg::Point3D& x ) {
   return std::sin( 4 * x[0] ) * std::sin( 8 * x[1] ) * std::sin( 2 * x[2] );
};
std::function< real_t( const hyteg::Point3D& ) > shellRhsU = []( const hyteg::Point3D& x ) {
   return 4 * std::sin( 8 * x[1] ) * std::sin( 2 * x[2] ) * std::cos( 4 * x[0] ) - 64 * std::cos( 4 * x[2] );
};
std::function< real_t( const hyteg::Point3D& ) > shellRhsV = []( const hyteg::Point3D& x ) {
   return 8 * std::sin( 4 * x[0] ) * std::sin( 2 * x[2] ) * std::cos( 8 * x[1] ) + 512 * std::cos( 8 * x[0] );
};
std::function< real_t( const hyteg::Point3D& ) > shellRhsW = []( const hyteg::Point3D& x ) {
   return 2 * std::sin( 4 * x[0] ) * std::sin( 8 * x[1] ) * std::cos( 2 * x[2] ) - 8 * std::cos( 2 * x[1] );
};

template < typename StokesOperator, typename FunctionType, typename VelocityMassOperator, bool ThreeDim >
void convergenceTest( uint_t level, real_t toleranceVelocityComponents, real_t tolerancePressure )
{
   const auto minLevel = 2;
   const auto maxLevel = level;

   const auto solverTargetResidual = 1e-8;
   const auto solverMaxIterations  = 10000;

   const auto enableVTK = false;

   WALBERLA_LOG_INFO_ON_ROOT( "Level: " << maxLevel );

   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();

   std::function< real_t( const Point3D& ) > exactU;
   std::function< real_t( const Point3D& ) > exactV;
   std::function< real_t( const Point3D& ) > exactW;
   std::function< real_t( const Point3D& ) > exactP;

   std::function< real_t( const Point3D& ) > rhsU;
   std::function< real_t( const Point3D& ) > rhsV;
   std::function< real_t( const Point3D& ) > rhsW;

   if ( ThreeDim )
   {
      meshInfo = MeshInfo::meshCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );

      exactU = shellExactU;
      exactV = shellExactV;
      exactW = shellExactW;
      exactP = shellExactP;

      rhsU = shellRhsU;
      rhsV = shellRhsV;
      rhsW = shellRhsW;
   }
   else
   {
      meshInfo = MeshInfo::meshRectangle( Point2D( { 0, 0 } ), Point2D( { 1, 1 } ), MeshInfo::CRISS, 1, 1 );

      exactU = plumeExactU;
      exactV = plumeExactV;
      exactP = plumeExactP;

      rhsU = plumeRhsU;
      rhsV = plumeRhsV;
   }

   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage );

   FunctionType r( "r", storage, minLevel, maxLevel );
   FunctionType f( "f", storage, minLevel, maxLevel );
   FunctionType u( "u", storage, minLevel, maxLevel );
   FunctionType uExact( "uExact", storage, minLevel, maxLevel );
   FunctionType err( "err", storage, minLevel, maxLevel );
   FunctionType tmp( "tmp", storage, minLevel, maxLevel );

   StokesOperator       A( storage, minLevel, maxLevel );
   VelocityMassOperator M( storage, minLevel, maxLevel );

   VTKOutput vtk( "../../output", "ElementwiseEpsilonMinResConvergenceTest", storage );
   vtk.add( u );
   vtk.add( uExact );
   vtk.add( f );
   vtk.add( err );

   if ( ThreeDim )
   {
     u.uvw().interpolate( { exactU, exactV, exactW }, maxLevel, hyteg::DirichletBoundary );
     uExact.uvw().interpolate( { exactU, exactV, exactW }, maxLevel );
     tmp.uvw().interpolate( { rhsU, rhsV, rhsW }, maxLevel );
   }
   else {
     u.uvw().interpolate( { exactU, exactV }, maxLevel, hyteg::DirichletBoundary );
     uExact.uvw().interpolate( { exactU, exactV }, maxLevel );
     tmp.uvw().interpolate( { rhsU, rhsV }, maxLevel );
   }

   uExact.p().interpolate( exactP, maxLevel );

   M.apply( tmp.uvw()[0], f.uvw()[0], maxLevel, All );
   M.apply( tmp.uvw()[1], f.uvw()[1], maxLevel, All );
   if ( ThreeDim )
   {
      M.apply( tmp.uvw()[2], f.uvw()[2], maxLevel, All );
   }

   std::shared_ptr< Solver< StokesOperator > > solver;

#ifdef HYTEG_BUILD_WITH_PETSC

   solver = std::make_shared< PETScBlockPreconditionedStokesSolver< StokesOperator > >(
       storage, maxLevel, solverTargetResidual, solverMaxIterations, 1, 1 );

#else

   solver = solvertemplates::stokesMinResSolver< StokesOperator >(
       storage, maxLevel, solverTargetResidual, solverMaxIterations, false );

#endif

   solver->solve( A, u, f, maxLevel );

   hyteg::vertexdof::projectMean( u.p(), maxLevel );
   hyteg::vertexdof::projectMean( uExact.p(), maxLevel );

   A.apply( u, r, maxLevel, hyteg::Inner | hyteg::NeumannBoundary );

   err.assign( { 1.0, -1.0 }, { u, uExact }, maxLevel );

   if ( enableVTK )
   {
      vtk.write( maxLevel );
   }

   uint_t velocityCompDoFs = numberOfGlobalDoFs< typename FunctionType::VelocityFunction_T::Tag >( *storage, maxLevel ) / 3;
   uint_t pressureCompDoFs = numberOfGlobalDoFs< typename FunctionType::PressureFunction_T::Tag >( *storage, maxLevel );

   real_t discr_l2_err_u = real_c(0);
   real_t discr_l2_err_v = real_c(0);
   real_t discr_l2_err_w = real_c(0);
   for( uint_t k = 0; k < err.uvw().getDimension(); k++ ) {
     discr_l2_err_u = std::sqrt( err.uvw()[k].dotGlobal( err.uvw()[k], maxLevel ) / real_c( velocityCompDoFs ) );
   }
   real_t discr_l2_err_p = std::sqrt( err.p().dotGlobal( err.p(), maxLevel ) / real_c( pressureCompDoFs ) );
   real_t residuum_l2    = std::sqrt( r.dotGlobal( r, maxLevel ) / real_c( 3 * velocityCompDoFs + pressureCompDoFs ) );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_u );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_v );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error w = " << discr_l2_err_w );
   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_p );
   WALBERLA_LOG_INFO_ON_ROOT( "residuum = " << residuum_l2 );

   WALBERLA_CHECK_LESS( discr_l2_err_u, toleranceVelocityComponents, "x component of velocity error is too large" )
   WALBERLA_CHECK_LESS( discr_l2_err_v, toleranceVelocityComponents, "y component of velocity error is too large" )
   WALBERLA_CHECK_LESS( discr_l2_err_w, toleranceVelocityComponents, "z component of velocity error is too large" )
   WALBERLA_CHECK_LESS( discr_l2_err_p, tolerancePressure, "pressure error is too large" )
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

#ifdef HYTEG_BUILD_WITH_PETSC
   hyteg::PETScManager petscManager( &argc, &argv );
#endif

   bool longrun = false;
   for ( int i = 0; i < argc; i++ )
   {
      auto arg = std::string( argv[i] );
      if ( arg == "--longrun" )
      {
         longrun = true;
      }
   }

   if ( longrun )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "--longrun flag was set!" );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "--longrun flag was NOT set!" );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" )
   WALBERLA_LOG_INFO_ON_ROOT( "### P1, 2D ###" )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   hyteg::convergenceTest< hyteg::P1P1ElementwiseAffineEpsilonStokesOperator,
                           hyteg::P1StokesFunction< walberla::real_t >,
                           hyteg::P1ElementwiseMassOperator,
                           false >( 3, 6.9e-2, 1.8 );
   hyteg::convergenceTest< hyteg::P1P1ElementwiseAffineEpsilonStokesOperator,
                           hyteg::P1StokesFunction< walberla::real_t >,
                           hyteg::P1ElementwiseMassOperator,
                           false >( 4, 2.7e-2, 1.7 );
   hyteg::convergenceTest< hyteg::P1P1ElementwiseAffineEpsilonStokesOperator,
                           hyteg::P1StokesFunction< walberla::real_t >,
                           hyteg::P1ElementwiseMassOperator,
                           false >( 5, 2.7e-2, 1.7 );

   WALBERLA_LOG_INFO_ON_ROOT( "" )
   WALBERLA_LOG_INFO_ON_ROOT( "### P1, 3D ###" )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   hyteg::convergenceTest< hyteg::P1P1ElementwiseAffineEpsilonStokesOperator,
                           hyteg::P1StokesFunction< walberla::real_t >,
                           hyteg::P1ElementwiseMassOperator,
                           true >( 3, 5.9e-1, 2.4e+1 );
   hyteg::convergenceTest< hyteg::P1P1ElementwiseAffineEpsilonStokesOperator,
                           hyteg::P1StokesFunction< walberla::real_t >,
                           hyteg::P1ElementwiseMassOperator,
                           true >( 4, 2.2e-1, 8.7e+0 );
   if ( longrun )
   {
      hyteg::convergenceTest< hyteg::P1P1ElementwiseAffineEpsilonStokesOperator,
                              hyteg::P1StokesFunction< walberla::real_t >,
                              hyteg::P1ElementwiseMassOperator,
                              true >( 5, 6.8e-2, 2.4e+0 );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" )
   WALBERLA_LOG_INFO_ON_ROOT( "### P2, 2D ###" )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   hyteg::convergenceTest< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator,
                           hyteg::P2P1TaylorHoodFunction< walberla::real_t >,
                           hyteg::P2ElementwiseMassOperator,
                           false >( 2, 2.7e-2, 2.2 );

   hyteg::convergenceTest< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator,
                           hyteg::P2P1TaylorHoodFunction< walberla::real_t >,
                           hyteg::P2ElementwiseMassOperator,
                           false >( 3, 1.8e-3, 3.8e-1 );

   hyteg::convergenceTest< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator,
                           hyteg::P2P1TaylorHoodFunction< walberla::real_t >,
                           hyteg::P2ElementwiseMassOperator,
                           false >( 4, 1.2e-4, 7.8e-2 );

   WALBERLA_LOG_INFO_ON_ROOT( "" )
   WALBERLA_LOG_INFO_ON_ROOT( "### P2, 3D ###" )
   WALBERLA_LOG_INFO_ON_ROOT( "" )

   hyteg::convergenceTest< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator,
                           hyteg::P2P1TaylorHoodFunction< walberla::real_t >,
                           hyteg::P2ElementwiseMassOperator,
                           true >( 2, 5.9e-1, 2.4e+1 );
   hyteg::convergenceTest< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator,
                           hyteg::P2P1TaylorHoodFunction< walberla::real_t >,
                           hyteg::P2ElementwiseMassOperator,
                           true >( 3, 2.2e-1, 8.7e+0 );

   if ( longrun )
   {
      hyteg::convergenceTest< hyteg::P2P1ElementwiseAffineEpsilonStokesOperator,
                              hyteg::P2P1TaylorHoodFunction< walberla::real_t >,
                              hyteg::P2ElementwiseMassOperator,
                              true >( 4, 6.8e-2, 2.4e+0 );
   }
}
