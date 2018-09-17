#include <tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp>
#include <tinyhhg_core/composites/P2P1TaylorHoodStokesOperator.hpp>
#include <tinyhhg_core/VTKWriter.hpp>
#include "core/mpi/MPIManager.h"

#include "tinyhhg_core/composites/P1StokesFunction.hpp"
#include "tinyhhg_core/composites/P1StokesOperator.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/DistributedBalancer.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"

#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/petsc/PETScLUSolver.hpp"
#include "tinyhhg_core/petsc/PETScVector.hpp"
#include "tinyhhg_core/petsc/PETScWrapper.hpp"

using walberla::real_t;

int main( int argc, char* argv[] )
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();

  std::string meshFileName = "../data/meshes/porous_fine.msh";

  hhg::MeshInfo              meshInfo = hhg::MeshInfo::fromGmshFile( meshFileName );
  hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  const uint_t level           = 2;
  const uint_t maxIterations   = 1000;
  const real_t targetResidual  = 1e-12;
  const bool   usePetsc        = true;

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage, timingTree );

#ifdef WALBERLA_BUILD_WITH_PARMETIS
  hhg::loadbalancing::distributed::parmetis( *storage );
#endif

  hhg::P2P1TaylorHoodFunction< real_t > r( "r", storage, level, level );
  hhg::P2P1TaylorHoodFunction< real_t > f( "f", storage, level, level );
  hhg::P2P1TaylorHoodFunction< real_t > u( "u", storage, level, level );

  hhg::P2P1TaylorHoodStokesOperator L( storage, level, level );

  std::function< real_t( const hhg::Point3D& ) > bc_x = []( const hhg::Point3D& x ) {
      if( x[0] < 1e-8 )
      {
        return 4.0 * x[1] * (1.0 - x[1]);
      } else
      {
        return 0.0;
      }
  };
  std::function< real_t( const hhg::Point3D& ) > rhs  = []( const hhg::Point3D& ) { return 0.0; };
  std::function< real_t( const hhg::Point3D& ) > zero = []( const hhg::Point3D& ) { return 0.0; };
  std::function< real_t( const hhg::Point3D& ) > ones = []( const hhg::Point3D& ) { return 1.0; };

  u.u.interpolate( bc_x, level, hhg::DirichletBoundary );
  u.v.interpolate( zero, level, hhg::DirichletBoundary );

  hhg::VTKOutput vtkOutput( "../output", "stokes_porous_taylor_hood" );

  vtkOutput.add( &r.u );
  vtkOutput.add( &r.v );
  vtkOutput.add( &r.p );

  vtkOutput.add( &f.u );
  vtkOutput.add( &f.v );
  vtkOutput.add( &f.p );

  vtkOutput.add( &u.u );
  vtkOutput.add( &u.v );
  vtkOutput.add( &u.p );

  timingTree->start( "Complete app" );

  vtkOutput.write( level, 0 );

#ifdef HHG_BUILD_WITH_PETSC
  if ( usePetsc )
  {
    PETScManager petscManager;
    f.u.interpolate(bc_x, level, hhg::DirichletBoundary);
    auto numerator = std::make_shared<hhg::P2P1TaylorHoodFunction<PetscInt>>("numerator", storage, level, level);
    numerator->enumerate(level);
    const uint_t localDoFs = hhg::numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
    const uint_t globalDoFs = hhg::numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, level );
    PETScLUSolver<real_t, hhg::P2P1TaylorHoodFunction, hhg::P2P1TaylorHoodStokesOperator> solver(numerator, localDoFs, globalDoFs);
    solver.solve( L, u, f, r, level, targetResidual, maxIterations, hhg::Inner | hhg::NeumannBoundary, true );
  }
  else
#else
  if ( usePetsc )
  {
    WALBERLA_LOG_INFO_ON_ROOT( "HHG was not built with PETSc - solving with default solver now..." );
  }
#endif
  {
    auto solver = hhg::MinResSolver< hhg::P2P1TaylorHoodFunction< real_t >,
                                     hhg::P2P1TaylorHoodStokesOperator      >( storage, level, level );

    solver.solve( L, u, f, r, level, targetResidual, maxIterations, hhg::Inner | hhg::NeumannBoundary, true );
  }

  vtkOutput.write( level, 1 );

  timingTree->stop( "Complete app" );
  WALBERLA_LOG_INFO_ON_ROOT( *timingTree );

  return EXIT_SUCCESS;
}
