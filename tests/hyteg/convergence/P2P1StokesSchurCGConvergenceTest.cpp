/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/StokesPCGSolverOld.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void P2P1SchurCGConvergenceTest( const uint_t & level, const MeshInfo & meshInfo )
{
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

  const uint_t minLevel = 2;

  hyteg::loadbalancing::roundRobin( setupStorage );

  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
  writeDomainPartitioningVTK( storage, "../../output", "P2P1Stokes2DSchurCGConvergence_Domain" );

  hyteg::P2P1TaylorHoodFunction< real_t >                      x( "x", storage, minLevel, level );
  hyteg::P2P1TaylorHoodFunction< real_t >                      x_exact( "x_exact", storage, minLevel, level );
  hyteg::P2P1TaylorHoodFunction< real_t >                      btmp( "btmp", storage, minLevel, level );
  hyteg::P2P1TaylorHoodFunction< real_t >                      b( "b", storage, minLevel, level );
  hyteg::P2P1TaylorHoodFunction< real_t >                      err( "err", storage, minLevel, level );
  hyteg::P2P1TaylorHoodFunction< real_t >                      residuum( "res", storage, minLevel, level );

  hyteg::P2P1TaylorHoodStokesOperator A( storage, minLevel, level );

  std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& xx ) { return real_c(20) * xx[0] * xx[1] * xx[1] * xx[1]; };
  std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& xx ) { return real_c(5) * xx[0] * xx[0] * xx[0] * xx[0] - real_c(5) * xx[1] * xx[1] * xx[1] * xx[1]; };
  std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& xx ) { return real_c(60) * std::pow( xx[0], 2.0 ) * xx[1] - real_c(20) * std::pow( xx[1], 3.0 ); };
  std::function< real_t( const hyteg::Point3D& ) > zero =   []( const hyteg::Point3D&    ) { return real_c(0); };

  x.uvw().interpolate( { exactU, exactV }, level, hyteg::DirichletBoundary );
  x_exact.uvw().interpolate({  exactU, exactV }, level );
  x_exact.p().interpolate( exactP, level );

//  VTKOutput vtkOutput("../../output", "P2P1Stokes2DSchurCGConvergence", storage);
//  vtkOutput.add( x.u );
//  vtkOutput.add( x.v );
//  vtkOutput.add( x.p() );
//  vtkOutput.add( x_exact.u );
//  vtkOutput.add( x_exact.v );
//  vtkOutput.add( x_exact.p() );
//  vtkOutput.add( err.u );
//  vtkOutput.add( err.v );
//  vtkOutput.add( err.p() );
//  vtkOutput.add( b.u );
//  vtkOutput.add( b.v );
//  vtkOutput.add( b.p() );
//  vtkOutput.write( level, 0 );

  uint_t localDoFs1 = hyteg::numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
  uint_t globalDoFs1 = hyteg::numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
  uint_t globalDoFsPerVelocityComponent = hyteg::numberOfGlobalDoFs< P2FunctionTag >( *storage, level );
  uint_t globalDoFsPressure = hyteg::numberOfGlobalDoFs< P1FunctionTag >( *storage, level );

  WALBERLA_LOG_INFO( "localDoFs1: " << localDoFs1 << " globalDoFs1: " << globalDoFs1 );

  auto coarseGrid     = std::make_shared< CGSolver< P2ConstantLaplaceOperator > >( storage, minLevel, minLevel, std::numeric_limits< uint_t >::max(), 1e-16 );
  auto smoother       = std::make_shared< GaussSeidelSmoother< P2ConstantLaplaceOperator > >();
  auto prolongation   = std::make_shared< P2toP2QuadraticProlongation >();
  auto restriction    = std::make_shared< P2toP2QuadraticRestriction >();
  auto velocitySolver = std::make_shared< GeometricMultigridSolver< P2ConstantLaplaceOperator > >( storage, smoother, coarseGrid, restriction, prolongation, minLevel, level, 2, 2 );
  auto loop           = std::make_shared< SolverLoop< P2ConstantLaplaceOperator > >( velocitySolver, 10 );
  StokesPCGSolverOld< P2P1TaylorHoodStokesOperator > solver( storage, loop, minLevel, level, 1e-10, 100, Inner | NeumannBoundary );

  walberla::WcTimer timer;
  solver.solve( A, x, b, level );
  timer.end();

  // hyteg::vertexdof::projectMean( x.p(), err.p(), level );
  // hyteg::vertexdof::projectMean( x_exact.p(), err.p(), level );

  WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
  A.apply( x, residuum, level, hyteg::Inner | hyteg::NeumannBoundary );

  err.assign( {1.0, -1.0}, {x, x_exact}, level );

  real_t discr_l2_err_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], level ) / (real_t) globalDoFsPerVelocityComponent );
  real_t discr_l2_err_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], level ) / (real_t) globalDoFsPerVelocityComponent );
  real_t discr_l2_err_p = std::sqrt( err.p().dotGlobal( err.p(), level ) / (real_t) globalDoFsPressure );
  real_t residuum_l2_1  = std::sqrt( residuum.dotGlobal( residuum, level, Inner ) / (real_t) globalDoFs1 );

  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_u );
  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_v );
  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_p );
  WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2_1 );

//  vtkOutput.write( level, 3 );

  WALBERLA_CHECK_LESS( discr_l2_err_u, 3.5e-03 );
  WALBERLA_CHECK_LESS( discr_l2_err_v, 2.4e-03 );
  WALBERLA_CHECK_LESS( discr_l2_err_p, 2.9e-01 );
  WALBERLA_CHECK_LESS( residuum_l2_1,  2.5e-09 );
}

}

using namespace hyteg;

int main( int argc, char* argv[] )
{
  walberla::Environment walberlaEnv( argc, argv );
  walberla::MPIManager::instance()->useWorldComm();

  P2P1SchurCGConvergenceTest( 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/quad_center_at_origin_4el.msh" ) );

  return EXIT_SUCCESS;
}

