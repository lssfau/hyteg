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

#include "hyteg/composites/P1P0StokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/misc/ExactStencilWeights.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScVersion.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/composites/P1P1StokesOperator.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

/// \param solverType
///     0: LU,
///     1: MinRes,
///     2: block preconditioned MinRes
void petscSolveTest( const uint_t & solverType, const uint_t & level, const MeshInfo & meshInfo, const real_t & resEps, const real_t & errEpsUSum, const real_t & errEpsP )
{
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

  hyteg::loadbalancing::roundRobin( setupStorage );

  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );
//  writeDomainPartitioningVTK( storage, "../../output", "P1P1StokesPetscSolve_Domain" );

  hyteg::P1P0StokesFunction< real_t >                      x( "x", storage, level, level );
  hyteg::P1P0StokesFunction< real_t >                      x_exact( "x_exact", storage, level, level );
  hyteg::P1P0StokesFunction< real_t >                      b( "b", storage, level, level );
  hyteg::P1P0StokesFunction< real_t >                      btmp( "btmp", storage, level, level );
  hyteg::P1P0StokesFunction< real_t >                      err( "err", storage, level, level );
  hyteg::P1P0StokesFunction< real_t >                      residuum( "res", storage, level, level );
  hyteg::P1P0StokesFunction< real_t >                      nullspace( "nullspace", storage, level, level );

  hyteg::P1P0StokesOperator A( storage, level, level, 0.1 );
  hyteg::P1ConstantMassOperator   M( storage, level, level );

  std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& xx ) { return - real_c(4) * std::cos( real_c(4) * xx[2] ); };
  std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& xx ) { return   real_c(8) * std::cos( real_c(8) * xx[0] ); };
  std::function< real_t( const hyteg::Point3D& ) > exactW = []( const hyteg::Point3D& xx ) { return - real_c(2) * std::cos( real_c(2) * xx[1] ); };

  std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& xx ) { return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ); };

  std::function< real_t( const hyteg::Point3D& ) > forceU = []( const hyteg::Point3D& xx ) { return 4*std::sin(8*xx[1])*std::sin(2*xx[2])*std::cos(4*xx[0])-64*std::cos(4*xx[2]); };
  std::function< real_t( const hyteg::Point3D& ) > forceV = []( const hyteg::Point3D& xx ) { return 8*std::sin(4*xx[0])*std::sin(2*xx[2])*std::cos(8*xx[1])+512*std::cos(8*xx[0]); };
  std::function< real_t( const hyteg::Point3D& ) > forceW = []( const hyteg::Point3D& xx ) { return 2*std::sin(4*xx[0])*std::sin(8*xx[1])*std::cos(2*xx[2])-8*std::cos(2*xx[1]); };

  std::function< real_t( const hyteg::Point3D& ) > zero =   []( const hyteg::Point3D&    ) { return real_c(0); };
  std::function< real_t( const hyteg::Point3D& ) > ones =   []( const hyteg::Point3D&    ) { return real_c(1); };


  btmp.uvw().interpolate( {forceU, forceV, forceW}, level, Inner );

  M.apply( btmp.uvw()[0], b.uvw()[0], level, All );
  M.apply( btmp.uvw()[1], b.uvw()[1], level, All );
  M.apply( btmp.uvw()[2], b.uvw()[2], level, All );


  b.uvw().interpolate( { exactU, exactV, exactW }, level, DirichletBoundary );
  b.p().interpolate( zero,   level, All );

  x.uvw().interpolate( { exactU, exactV, exactW }, level, DirichletBoundary );

  x_exact.uvw().interpolate( { exactU, exactV, exactW }, level );
  x_exact.p().interpolate( exactP, level );

  nullspace.p().interpolate( ones, level, All );

//  VTKOutput vtkOutput("../../output", "P1P1StokesPetscSolve", storage);
//  vtkOutput.add( x.u );
//  vtkOutput.add( x.v );
//  vtkOutput.add( x.w );
//  vtkOutput.add( x.p );
//  vtkOutput.add( x_exact.u );
//  vtkOutput.add( x_exact.v );
//  vtkOutput.add( x_exact.w );
//  vtkOutput.add( x_exact.p );
//  vtkOutput.add( err.u );
//  vtkOutput.add( err.v );
//  vtkOutput.add( err.w );
//  vtkOutput.add( err.p );
//  vtkOutput.add( b.u );
//  vtkOutput.add( b.v );
//  vtkOutput.add( b.w );
//  vtkOutput.add( b.p );
//  vtkOutput.write( level, 0 );

  uint_t localDoFs1 = hyteg::numberOfLocalDoFs< P1StokesFunctionTag >( *storage, level );
  uint_t globalDoFs1 = hyteg::numberOfGlobalDoFs< P1StokesFunctionTag >( *storage, level );

  WALBERLA_LOG_INFO( "localDoFs: " << localDoFs1 << " globalDoFs: " << globalDoFs1 );

  PETScMinResSolver< P1P0StokesOperator > solver_1( storage, level );

  walberla::WcTimer timer;

    WALBERLA_LOG_INFO_ON_ROOT( "MinRes solver ..." )
      solver_1.solve( A, x, b, level );
      // The second solve() call is intended!
      // The solvers / assembly procedures must be tested for repeated solves
      // (there have already been numerous issues).



  timer.end();

  hyteg::dg::projectMean( x.p(), level );
  hyteg::dg::projectMean( x_exact.p(), level );

  WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
  A.apply( x, residuum, level, hyteg::Inner );

  err.assign( {1.0, -1.0}, {x, x_exact}, level );

  real_t discr_l2_err_1_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], level ) / (real_t) globalDoFs1 );
  real_t discr_l2_err_1_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], level ) / (real_t) globalDoFs1 );
  real_t discr_l2_err_1_w = std::sqrt( err.uvw()[2].dotGlobal( err.uvw()[2], level ) / (real_t) globalDoFs1 );
  real_t discr_l2_err_1_p = std::sqrt( err.p().dotGlobal( err.p(), level ) / (real_t) globalDoFs1 );
  real_t residuum_l2_1  = std::sqrt( residuum.dotGlobal( residuum, level ) / (real_t) globalDoFs1 );

  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_1_u );
  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_1_v );
  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error w = " << discr_l2_err_1_w );
  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_1_p );
  WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2_1 );

//  vtkOutput.write( level, 1 );

/*  WALBERLA_CHECK_LESS( residuum_l2_1, resEps );
  WALBERLA_CHECK_LESS( discr_l2_err_1_u + discr_l2_err_1_v + discr_l2_err_1_w, errEpsUSum );
  WALBERLA_CHECK_LESS( discr_l2_err_1_p, errEpsP);*/
}

} 

using namespace hyteg;

int main( int argc, char* argv[] )
{
    walberla::MPIManager::instance()->initializeMPI(&argc, &argv);
    walberla::MPIManager::instance()->useWorldComm();
    hyteg::PETScManager petscManager(&argc, &argv);


  petscSolveTest( 1, 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_center_at_origin_24el.msh" ), 8.0e-15, 0.118, 2.78653 );
    petscSolveTest( 1, 4, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_center_at_origin_24el.msh" ), 8.0e-15, 0.118, 2.78653 );
    petscSolveTest( 1, 5, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_center_at_origin_24el.msh" ), 8.0e-15, 0.118, 2.78653 );

  return EXIT_SUCCESS;
}
