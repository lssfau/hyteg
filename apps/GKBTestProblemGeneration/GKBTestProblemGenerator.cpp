/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include "hyteg/communication/Syncing.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/petsc/PETScExportFunctionAsVector.hpp"
#include "hyteg/petsc/PETScExportLinearSystem.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

using walberla::math::pi;

std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& x ) {
   return std::sin( 2 * pi * x[0] ) * std::cos( pi * x[1] );
};

std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& x ) {
   return -2.0 * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
};

std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& x ) {
   return 2.5 * pi * std::cos( 2 * pi * x[0] ) * std::cos( pi * x[1] );
};
std::function< real_t( const hyteg::Point3D& ) > rhsU = []( const hyteg::Point3D& ) { return 0; };
std::function< real_t( const hyteg::Point3D& ) > rhsV = []( const hyteg::Point3D& x ) {
   return -12.5 * pi * pi * std::cos( 2 * pi * x[0] ) * std::sin( pi * x[1] );
};

void runBenchmark( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petscManager( &argc, &argv );

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./GKBTestProblemGenerator.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }
   /////////////// Parameters ///////////////
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );
   const uint_t                        level    = mainConf.getParameter< uint_t >( "level" );
   const bool                          writeVTK = mainConf.getParameter< bool >( "vtk" );
   ;

   auto meshInfo = MeshInfo::meshRectangle( Point2D( {0, 0} ), Point2D( {1, 1} ), MeshInfo::CRISS, 1, 1 );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   writeDomainPartitioningVTK( storage, "./output", "GKBTestProblemDomain" );

   P2P1TaylorHoodFunction< real_t > r( "r", storage, level, level );
   P2P1TaylorHoodFunction< real_t > f( "f", storage, level, level );
   P2P1TaylorHoodFunction< real_t > u( "u", storage, level, level );
   P2P1TaylorHoodFunction< real_t > Au( "Au", storage, level, level );
   P2P1TaylorHoodFunction< real_t > u_exact( "u_exact", storage, level, level );
   P2P1TaylorHoodFunction< real_t > err( "err", storage, level, level );

   typedef P2P1TaylorHoodStokesOperator StokesOperator;
   typedef P2ConstantMassOperator       MassOperator;

   StokesOperator L( storage, level, level );
   MassOperator   M( storage, level, level );

   u.uvw().interpolate( { exactU, exactV }, level, DirichletBoundary );

   Au.uvw().interpolate( { rhsU, rhsV }, level, All );

   M.apply( Au.uvw()[0], f.uvw()[0], level, All );
   M.apply( Au.uvw()[1], f.uvw()[1], level, All );

   Au.uvw()[0].setToZero( level );
   Au.uvw()[1].setToZero( level );
   Au.p().setToZero( level );

   u_exact.uvw().interpolate( { exactU, exactV }, level, All );
   u_exact.p().interpolate( exactP, level, All );

   vertexdof::projectMean( u_exact.p(), level );

   communication::syncVectorFunctionBetweenPrimitives( u_exact.uvw(), level );
   communication::syncFunctionBetweenPrimitives( u_exact.p(), level );

   auto solver = std::make_shared< PETScLUSolver< StokesOperator > >( storage, level );

   const uint_t numPoints         = numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
   const uint_t numPointsVelocity = numberOfGlobalDoFs< P2FunctionTag >( *storage, level );
   const uint_t numPointsPressure = numberOfGlobalDoFs< P1FunctionTag >( *storage, level );

   WALBERLA_LOG_INFO_ON_ROOT( "level:                                    " << level )

   WALBERLA_LOG_INFO_ON_ROOT( "num unknowns total:                       " << numPoints )
   WALBERLA_LOG_INFO_ON_ROOT( "num velocity unknowns (u and v combined): " << 2 * numPointsVelocity )
   WALBERLA_LOG_INFO_ON_ROOT( "num pressure unknowns:                    " << numPointsPressure )
   real_t currRes = 0, oldRes = 0;

   L.apply( u, Au, level, Inner | NeumannBoundary );
   r.assign( {1.0, -1.0}, {f, Au}, level, Inner | NeumannBoundary );
   oldRes = std::sqrt( r.dotGlobal( r, level, All ) ) / real_c( numPoints );

   err.assign( {1.0, -1.0}, {u, u_exact}, level, All );

   auto discr_l2_err_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], level, All ) / real_c( numPointsVelocity ) );
   auto discr_l2_err_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], level, All ) / real_c( numPointsVelocity ) );
   auto discr_l2_err_p = std::sqrt( err.p().dotGlobal( err.p(), level, All ) / real_c( numPointsPressure ) );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
       "%15s|%15s|%15s|%15s|%15s|%15s", "iteration", "residual", "residual redct", "L2 error u", "L2 error v", "L2 error p" ) )
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
       "%15s|%15e|%15s|%15e|%15e|%15e", "initial", oldRes, "-", discr_l2_err_u, discr_l2_err_v, discr_l2_err_p ) )

   VTKOutput vtkOutput( "./output", "GKBTestProblem", storage );
   vtkOutput.add( u );
   vtkOutput.add( u_exact );
   vtkOutput.add( err );

   if ( writeVTK )
   {
      vtkOutput.write( level, 0 );
   }

   solver->solve( L, u, f, level );

   vertexdof::projectMean( u.p(), level );

   L.apply( u, Au, level, Inner | NeumannBoundary );
   r.assign( {1.0, -1.0}, {f, Au}, level, Inner | NeumannBoundary );
   currRes = std::sqrt( r.dotGlobal( r, level, All ) ) / real_c( numPoints );

   auto residualReduction = currRes / oldRes;
   oldRes                 = currRes;

   err.assign( {1.0, -1.0}, {u, u_exact}, level );
   discr_l2_err_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], level, All ) / real_c( numPointsVelocity ) );
   discr_l2_err_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], level, All ) / real_c( numPointsVelocity ) );
   discr_l2_err_p = std::sqrt( err.p().dotGlobal( err.p(), level, All ) / real_c( numPointsPressure ) );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format(
       "%15s|%15e|%15e|%15e|%15e|%15e", "solved", currRes, residualReduction, discr_l2_err_u, discr_l2_err_v, discr_l2_err_p ) )

   if ( writeVTK )
   {
      vtkOutput.write( level, 1 );
   }

   exportLinearSystem< P2P1TaylorHoodStokesOperator, P2P1TaylorHoodFunction, P2P1TaylorHoodFunctionTag >(
       L,
       f,
       u_exact,
       "./output/stokes_matrix_level_" + std::to_string( level ) + ".m",
       "L_level_" + std::to_string( level ),
       "./output/stokes_f_level_" + std::to_string( level ) + ".m",
       "f_level_" + std::to_string( level ),
       storage,
       level,
       true,
       false );

   exportFunction< P2P1TaylorHoodFunction, P2P1TaylorHoodFunctionTag >( u,
                                                                        "./output/stokes_u_level_" + std::to_string( level ) +
                                                                            ".m",
                                                                        "u_level_" + std::to_string( level ),
                                                                        storage,
                                                                        level,
                                                                        false );

   exportFunction< P2P1TaylorHoodFunction, P2P1TaylorHoodFunctionTag >( u_exact,
                                                                        "./output/stokes_exact_level_" + std::to_string( level ) +
                                                                            ".m",
                                                                        "exact_level_" + std::to_string( level ),
                                                                        storage,
                                                                        level,
                                                                        false );
}
int main( int argc, char* argv[] )
{
   runBenchmark( argc, argv );

   return 0;
}
