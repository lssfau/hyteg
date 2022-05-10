#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/misc/ExactStencilWeights.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/solvers/preconditioners/stokes/StokesPressureBlockPreconditioner.hpp"
#include "hyteg/solvers/MinresSolver.hpp"
#include "hyteg/elementwiseoperators/P2P1ElementwiseAffineEpsilonStokesOperator.hpp"
#include "hyteg/petsc/PETScBlockPreconditionedStokesSolver.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

// SolVi benchmark (Circular inclusion) 
/*  Described in section 5.3 of [1] with analytical solution from [2] (equations 26, 34).
The viscosity contains a jump from visc_matrix to visc_inclusion in a circle located at the center of the domain. 
Difficulty for FE methods: the mesh can not align with the viscosity jump.

[1]: "On the choice of finite element for applications in geodynamics" 
by Cedric Thieulot and Wolfgang Bangerth

[2]: "Analytical solutions for deformable elliptical inclusions in general shear" 
by Daniel W. Schmid and Yuri Yu. Podladchikov
*/

void SolViBenchmark( const uint_t & level, const uint_t & nxy, const real_t r_inclusion, const real_t & visc_inclusion, const real_t & visc_matrix)
{
  // storage and domain
  auto meshInfo = MeshInfo::meshRectangle( Point2D( {-1, -1} ), Point2D( {1, 1} ), MeshInfo::CRISSCROSS, nxy, nxy );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );  
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  hyteg::loadbalancing::roundRobin( setupStorage );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );
  writeDomainPartitioningVTK( storage, "../../output", "SolViBenchmark_Domain" );


  // function setup
  hyteg::P2P1TaylorHoodFunction< real_t >                      x( "x", storage, level, level );
  hyteg::P2P1TaylorHoodFunction< real_t >                      x_exact( "x_exact", storage, level, level );
  hyteg::P2P1TaylorHoodFunction< real_t >                      btmp( "btmp", storage, level, level );
  hyteg::P2P1TaylorHoodFunction< real_t >                      b( "b", storage, level, level );
  hyteg::P2P1TaylorHoodFunction< real_t >                      err( "err", storage, level, level );
  hyteg::P2P1TaylorHoodFunction< real_t >                      residuum( "res", storage, level, level );

  // radius helper function
  std::function< real_t( const hyteg::Point3D& ) > rad = []( const hyteg::Point3D& xx ) { return sqrt(std::pow(xx[0],2.0) + std::pow(xx[1],2.0)); };
 
  // viscosity function and operator setup
  std::function< real_t( const hyteg::Point3D& ) > viscosity = [r_inclusion, visc_inclusion, visc_matrix, rad](const hyteg::Point3D& xx ) { 
    if(rad(xx) <= r_inclusion)  
        return visc_inclusion;
    else return visc_matrix;
  };
  hyteg::P2P1ElementwiseAffineEpsilonStokesOperator A( storage, level, level, viscosity);


  // exact solution
  const real_t C_visc = visc_matrix/(visc_inclusion + visc_matrix);
  std::function< real_t( const hyteg::Point3D& ) > exactU = [C_visc, r_inclusion]( const hyteg::Point3D& xx ) { return real_c(2) * C_visc * xx[0]; };
  std::function< real_t( const hyteg::Point3D& ) > exactV = [C_visc, r_inclusion]( const hyteg::Point3D& xx ) { return real_c(-2) * C_visc * xx[1]; };
  std::function< real_t( const hyteg::Point3D& ) > exactP = [C_visc, r_inclusion, visc_inclusion, visc_matrix, rad]( const hyteg::Point3D& xx ) { 
    if(rad(xx) <= r_inclusion) return real_c(0);
    else return real_c(4) * C_visc * (visc_inclusion - visc_matrix) * std::pow(r_inclusion, 2.0)/(std::pow(xx[0],2.0) + std::pow(xx[1],2.0)) * cos(2*atan2(xx[1],xx[0])); 
  };
  x_exact.uvw().interpolate( { exactU, exactV }, level );
  x_exact.p().interpolate( exactP, level );
  x.uvw().interpolate( { exactU, exactV }, level, hyteg::DirichletBoundary );
  x.p().interpolate( { exactP }, level, hyteg::DirichletBoundary );
   

  // Right-hand-side
  A.apply( x_exact, b, level, hyteg::Inner | hyteg::NeumannBoundary );

 
  // Visualization
  VTKOutput vtkOutput("../../output", "SolViBenchmark", storage);
  vtkOutput.add( x.uvw() );
  vtkOutput.add( x.p() );
  vtkOutput.add( x_exact.uvw() );
  vtkOutput.add( x_exact.p() );
  //vtkOutput.add( err.uvw() );
  //vtkOutput.add( err.p() );
  //vtkOutput.add( b.uvw() );
  //vtkOutput.add( b.p() );
  
  vtkOutput.write( level, 0 );
  uint_t localDoFs1 = hyteg::numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
  uint_t globalDoFs1 = hyteg::numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
  WALBERLA_LOG_INFO( "localDoFs1: " << localDoFs1 << " globalDoFs1: " << globalDoFs1 );
  

  // Initial errors and residual
  A.apply( x, btmp, level, hyteg::Inner | hyteg::NeumannBoundary );
  residuum.assign({1.0, -1.0}, {b, btmp}, level, hyteg::Inner | hyteg::NeumannBoundary);
  err.assign( {1.0, -1.0}, {x, x_exact}, level, hyteg::Inner | hyteg::NeumannBoundary );
  real_t discr_l2_err_1_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], level, hyteg::Inner | hyteg::NeumannBoundary  ) / (real_t) globalDoFs1 );
  real_t discr_l2_err_1_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], level, hyteg::Inner | hyteg::NeumannBoundary ) / (real_t) globalDoFs1 );
  real_t discr_l2_err_1_p = std::sqrt( err.p().dotGlobal( err.p(), level, hyteg::Inner | hyteg::NeumannBoundary ) / (real_t) globalDoFs1 );
  real_t residuum_l2_1  = std::sqrt( residuum.dotGlobal( residuum, level, hyteg::Inner | hyteg::NeumannBoundary ) / (real_t) globalDoFs1 );
  WALBERLA_LOG_INFO_ON_ROOT( "initial errors and residual:");
  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_1_u );
  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_1_v );
  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_1_p );
  WALBERLA_LOG_INFO_ON_ROOT( "residual = " << residuum_l2_1 );
  
  // Solve
  //auto solver = hyteg::MinResSolver< P2P1ElementwiseAffineEpsilonStokesOperator >( storage, level, level );
  PETScLUSolver< P2P1ElementwiseAffineEpsilonStokesOperator > solver( storage, level );
  /* auto solver = hyteg::PETScBlockPreconditionedStokesSolver< P2P1ElementwiseAffineEpsilonStokesOperator >( 
    storage,
    level,
    1e-12,
    std::numeric_limits< PetscInt >::max(),
    3,
    1
  );*/
  //solver.setPrintInfo(true);
  walberla::WcTimer timer;
  solver.solve( A, x, b, level );
  timer.end();

  hyteg::vertexdof::projectMean( x.p(), level );
  hyteg::vertexdof::projectMean( x_exact.p(), level );


  WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
  A.apply( x, btmp, level, hyteg::Inner | hyteg::NeumannBoundary );
  residuum.assign({1.0, -1.0}, {b, btmp}, level);
  err.assign( {1.0, -1.0}, {x, x_exact}, level );
  discr_l2_err_1_u = std::sqrt( err.uvw()[0].dotGlobal( err.uvw()[0], level, hyteg::Inner | hyteg::NeumannBoundary  ) / (real_t) globalDoFs1 );
  discr_l2_err_1_v = std::sqrt( err.uvw()[1].dotGlobal( err.uvw()[1], level, hyteg::Inner | hyteg::NeumannBoundary ) / (real_t) globalDoFs1 );
  discr_l2_err_1_p = std::sqrt( err.p().dotGlobal( err.p(), level, hyteg::Inner | hyteg::NeumannBoundary ) / (real_t) globalDoFs1 );
  residuum_l2_1  = std::sqrt( residuum.dotGlobal( residuum, level, hyteg::Inner | hyteg::NeumannBoundary ) / (real_t) globalDoFs1 );
   WALBERLA_LOG_INFO_ON_ROOT( "final errors and residual:");
  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error u = " << discr_l2_err_1_u );
  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error v = " << discr_l2_err_1_v );
  WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error p = " << discr_l2_err_1_p );
  WALBERLA_LOG_INFO_ON_ROOT( "residual = " << residuum_l2_1 );

  vtkOutput.write( level, 1 );
}

}

using namespace hyteg;

int main( int argc, char* argv[] )
{
  walberla::Environment walberlaEnv( argc, argv );
  walberla::MPIManager::instance()->useWorldComm();
  PETScManager petscManager( &argc, &argv );

  SolViBenchmark( 3, 10, 0.2, 1000, 1);

  return EXIT_SUCCESS;
}

