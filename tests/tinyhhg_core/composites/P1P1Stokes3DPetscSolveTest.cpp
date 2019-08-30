#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/misc/ExactStencilWeights.hpp"
#include "tinyhhg_core/p1functionspace/P1ConstantOperator.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/petsc/PETScLUSolver.hpp"
#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/composites/P1StokesOperator.hpp"

#ifndef HHG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHHG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void petscSolveTest( const uint_t & level, const MeshInfo & meshInfo, const real_t & resEps, const real_t & errEpsUSum, const real_t & errEpsP )
{
  PETScManager petscManager;

  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

  hyteg::loadbalancing::roundRobin( setupStorage );

  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
  writeDomainPartitioningVTK( storage, "../../output", "P1P1StokesPetscSolve_Domain" );

  hyteg::P1StokesFunction< real_t >                      x( "x", storage, level, level );
  hyteg::P1StokesFunction< real_t >                      x_exact( "x_exact", storage, level, level );
  hyteg::P1StokesFunction< real_t >                      b( "b", storage, level, level );
  hyteg::P1StokesFunction< real_t >                      btmp( "btmp", storage, level, level );
  hyteg::P1StokesFunction< real_t >                      err( "err", storage, level, level );
  hyteg::P1StokesFunction< real_t >                      residuum( "res", storage, level, level );
  hyteg::P1StokesFunction< real_t >                      nullspace( "nullspace", storage, level, level );

  hyteg::P1StokesOperator A( storage, level, level );
  hyteg::P1ConstantMassOperator   M( storage, level, level );

#if 0
  std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& xx ) { return - real_c(4) * std::cos( real_c(4) * xx[2] ); };
  std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& xx ) { return   real_c(8) * std::cos( real_c(8) * xx[0] ); };
  std::function< real_t( const hyteg::Point3D& ) > exactW = []( const hyteg::Point3D& xx ) { return - real_c(2) * std::cos( real_c(2) * xx[1] ); };

  std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& xx ) { return std::sin( 4 * xx[0] ) * std::sin( 8 * xx[1] ) * std::sin( 2 * xx[2] ); };

  std::function< real_t( const hyteg::Point3D& ) > forceU = []( const hyteg::Point3D& xx ) { return 4*std::sin(8*xx[1])*std::sin(2*xx[2])*std::cos(4*xx[0])-64*std::cos(4*xx[2]); };
  std::function< real_t( const hyteg::Point3D& ) > forceV = []( const hyteg::Point3D& xx ) { return 8*std::sin(4*xx[0])*std::sin(2*xx[2])*std::cos(8*xx[1])+512*std::cos(8*xx[0]); };
  std::function< real_t( const hyteg::Point3D& ) > forceW = []( const hyteg::Point3D& xx ) { return 2*std::sin(4*xx[0])*std::sin(8*xx[1])*std::cos(2*xx[2])-8*std::cos(2*xx[1]); };

  std::function< real_t( const hyteg::Point3D& ) > zero =   []( const hyteg::Point3D&    ) { return real_c(0); };
  std::function< real_t( const hyteg::Point3D& ) > ones =   []( const hyteg::Point3D&    ) { return real_c(1); };

#endif

  std::function< real_t( const hyteg::Point3D& ) > exactU = []( const hyteg::Point3D& xx ) { return real_c(20) * xx[0] * xx[1] * xx[1] * xx[1]; };
  std::function< real_t( const hyteg::Point3D& ) > exactV = []( const hyteg::Point3D& xx ) { return real_c(5) * xx[0] * xx[0] * xx[0] * xx[0] - real_c(5) * xx[1] * xx[1] * xx[1] * xx[1]; };
  std::function< real_t( const hyteg::Point3D& ) > exactW = []( const hyteg::Point3D&    ) { return real_c(0); };
  std::function< real_t( const hyteg::Point3D& ) > exactP = []( const hyteg::Point3D& xx ) { return real_c(60) * std::pow( xx[0], 2.0 ) * xx[1] - real_c(20) * std::pow( xx[1], 3.0 ); };
  std::function< real_t( const hyteg::Point3D& ) > zero =   []( const hyteg::Point3D&    ) { return real_c(0); };
  std::function< real_t( const hyteg::Point3D& ) > ones =   []( const hyteg::Point3D&    ) { return real_c(1); };

#if 0
  btmp.u.interpolate( forceU, level, Inner );
  btmp.v.interpolate( forceV, level, Inner );
  btmp.w.interpolate( forceW, level, Inner );

  M.apply( btmp.u, b.u, level, All );
  M.apply( btmp.v, b.v, level, All );
  M.apply( btmp.w, b.w, level, All );

#endif

  b.u.interpolate( exactU, level, DirichletBoundary );
  b.v.interpolate( exactV, level, DirichletBoundary );
  b.w.interpolate( exactW, level, DirichletBoundary );
  b.p.interpolate( zero,   level, All );

  x.u.interpolate( exactU, level, DirichletBoundary );
  x.v.interpolate( exactV, level, DirichletBoundary );
  x.w.interpolate( exactW, level, DirichletBoundary );

  x_exact.u.interpolate( exactU, level );
  x_exact.v.interpolate( exactV, level );
  x_exact.w.interpolate( exactW, level );
  x_exact.p.interpolate( exactP, level );

  nullspace.p.interpolate( ones, level, All );

  VTKOutput vtkOutput("../../output", "P1P1StokesPetscSolve", storage);
  vtkOutput.add( x.u );
  vtkOutput.add( x.v );
  vtkOutput.add( x.w );
  vtkOutput.add( x.p );
  vtkOutput.add( x_exact.u );
  vtkOutput.add( x_exact.v );
  vtkOutput.add( x_exact.w );
  vtkOutput.add( x_exact.p );
  vtkOutput.add( err.u );
  vtkOutput.add( err.v );
  vtkOutput.add( err.w );
  vtkOutput.add( err.p );
  vtkOutput.add( b.u );
  vtkOutput.add( b.v );
  vtkOutput.add( b.w );
  vtkOutput.add( b.p );
  vtkOutput.write( level, 0 );

  uint_t localDoFs1 = hyteg::numberOfLocalDoFs< P1StokesFunctionTag >( *storage, level );
  uint_t globalDoFs1 = hyteg::numberOfGlobalDoFs< P1StokesFunctionTag >( *storage, level );

  WALBERLA_LOG_INFO( "localDoFs1: " << localDoFs1 << " globalDoFs1: " << globalDoFs1 );

  PETScLUSolver< hyteg::P1StokesOperator > solver_1( storage, level );

  walberla::WcTimer timer;
  solver_1.solve( A, x, b, level );
  timer.end();

  hyteg::vertexdof::projectMean( x.p, level );
  hyteg::vertexdof::projectMean( x_exact.p, level );

  WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
  A.apply( x, residuum, level, hyteg::Inner );

  err.assign( {1.0, -1.0}, {x, x_exact}, level );

  real_t discr_l2_err_1_u = std::sqrt( err.u.dotGlobal( err.u, level ) / (real_t) globalDoFs1 );
  real_t discr_l2_err_1_v = std::sqrt( err.v.dotGlobal( err.v, level ) / (real_t) globalDoFs1 );
  real_t discr_l2_err_1_w = std::sqrt( err.w.dotGlobal( err.w, level ) / (real_t) globalDoFs1 );
  real_t discr_l2_err_1_p = std::sqrt( err.p.dotGlobal( err.p, level ) / (real_t) globalDoFs1 );
  real_t residuum_l2_1  = std::sqrt( residuum.dotGlobal( residuum, level ) / (real_t) globalDoFs1 );

  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_1_u );
  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_1_v );
  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error w = " << discr_l2_err_1_w );
  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_1_p );
  WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2_1 );

  vtkOutput.write( level, 1 );

  WALBERLA_CHECK_LESS( residuum_l2_1, resEps );
  WALBERLA_CHECK_LESS( discr_l2_err_1_u + discr_l2_err_1_v + discr_l2_err_1_w, errEpsUSum );
  WALBERLA_CHECK_LESS( discr_l2_err_1_p, errEpsP);
}

}

using namespace hyteg;

int main( int argc, char* argv[] )
{
  walberla::Environment walberlaEnv( argc, argv );
  walberla::MPIManager::instance()->useWorldComm();

  petscSolveTest( 3, hyteg::MeshInfo::fromGmshFile( "../../data/meshes/3D/cube_center_at_origin_24el.msh" ), 8.0e-15, 0.118, 2.78653 );

  return EXIT_SUCCESS;
}
