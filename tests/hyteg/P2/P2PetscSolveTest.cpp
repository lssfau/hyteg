#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/VTKWriter.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/misc/ExactStencilWeights.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/FunctionProperties.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
WALBERLA_ABORT( "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON" )
#endif

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

void petscSolveTest( const uint_t & level, const std::string & meshFileName, const real_t & errEps )
{
   WALBERLA_LOG_INFO_ON_ROOT( "##### Mesh file: " << meshFileName << " / level: " << level << " #####" )

   PETScManager petscManager;

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFileName );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   hyteg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   hyteg::P2Function< real_t >                      x( "x", storage, level, level + 1 );
   hyteg::P2Function< real_t >                      x_exact( "x_exact", storage, level, level + 1 );
   hyteg::P2Function< real_t >                      b( "x", storage, level, level + 1 );
   hyteg::P2Function< real_t >                      err( "err", storage, level, level + 1 );
   hyteg::P2Function< real_t >                      residuum( "err", storage, level, level + 1 );

   hyteg::P2ConstantLaplaceOperator A( storage, level, level + 1 );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) { return sin( xx[0] ) * sinh( xx[1] ); };
   walberla::math::seedRandomGenerator( 0 );
   std::function< real_t( const Point3D& ) > rand = []( const Point3D& ) { return walberla::math::realRandom( 0.0, 1.0 ); };

   x.interpolate( exact, level, hyteg::DirichletBoundary );
   x.interpolate( rand, level, hyteg::Inner );
   b.interpolate( exact, level, hyteg::DirichletBoundary );
   x_exact.interpolate( exact, level );

//   x.interpolate( exact, level + 1, hyteg::DirichletBoundary );
//   x.interpolate( rand, level + 1, hyteg::Inner );
//   b.interpolate( exact, level + 1, hyteg::DirichletBoundary );
//   x_exact.interpolate( exact, level + 1 );

   uint_t localDoFs1 = hyteg::numberOfLocalDoFs< P2FunctionTag >( *storage, level );
//   uint_t localDoFs2 = hyteg::numberOfLocalDoFs< P2FunctionTag >( *storage, level + 1 );
   uint_t globalDoFs1 = hyteg::numberOfGlobalDoFs< P2FunctionTag >( *storage, level );
//   uint_t globalDoFs2 = hyteg::numberOfGlobalDoFs< P2FunctionTag >( *storage, level + 1 );

   WALBERLA_LOG_INFO( "localDoFs1: " << localDoFs1 << " globalDoFs1: " << globalDoFs1 );
//   WALBERLA_LOG_INFO( "localDoFs2: " << localDoFs2 << " globalDoFs2: " << globalDoFs2 );

   PETScLUSolver< hyteg::P2ConstantLaplaceOperator > solver_1( storage, level );
//   PETScLUSolver< real_t, hyteg::P2Function, hyteg::P2ConstantLaplaceOperator > solver_2( numerator, localDoFs2, globalDoFs2 );

   walberla::WcTimer timer;
   solver_1.solve( A, x, b, level );
//   solver_2.solve( A, x, b, x, level + 1, 0, 0 );
   timer.end();

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   A.apply( x, residuum, level, hyteg::Inner );
//   A.apply( x, residuum, level + 1, hyteg::Inner );

   err.assign( {1.0, -1.0}, {x, x_exact}, level );
//   err.assign( {1.0, -1.0}, {x, x_exact}, level + 1 );

   real_t discr_l2_err_1 = std::sqrt( err.dotGlobal( err, level ) / (real_t) globalDoFs1 );
//   real_t discr_l2_err_2 = std::sqrt( err.dotGlobal( err, level + 1 ) / (real_t) globalDoFs2 );
   real_t residuum_l2_1  = std::sqrt( residuum.dotGlobal( residuum, level ) / (real_t) globalDoFs1 );
//   real_t residuum_l2_2  = std::sqrt( residuum.dotGlobal( residuum, level + 1 ) / (real_t) globalDoFs2 );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error 1 = " << discr_l2_err_1 );
//   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error 2 = " << discr_l2_err_2 );
//   WALBERLA_LOG_INFO_ON_ROOT( "error ratio = " << ( discr_l2_err_1 / discr_l2_err_2 ) );
   WALBERLA_LOG_INFO_ON_ROOT( "residuum 1 = " << residuum_l2_1 );
//   WALBERLA_LOG_INFO_ON_ROOT( "residuum 2 = " << residuum_l2_2 );

   VTKOutput vtkOutput("../../output", "P2PetscSolve", storage);
   vtkOutput.add( x );
   vtkOutput.add( x_exact );
   vtkOutput.add( err );
   vtkOutput.add( residuum );
   vtkOutput.write( level );

   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( residuum_l2_1, 0.0, 1e-15 );
   //WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( residuum_l2_2, 0.0, 1e-15 );
   //WALBERLA_CHECK_LESS( 8.0, ( discr_l2_err_1 / discr_l2_err_2 ) );

   WALBERLA_CHECK_LESS( discr_l2_err_1, errEps );

}

}

using namespace hyteg;

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   petscSolveTest( 3, "../../data/meshes/quad_4el.msh",                  3.0e-07 );
   petscSolveTest( 3, "../../data/meshes/3D/tet_1el.msh",                3.0e-07 );
   petscSolveTest( 3, "../../data/meshes/3D/pyramid_2el.msh",            2.7e-06 );
   petscSolveTest( 3, "../../data/meshes/3D/pyramid_4el.msh",            3.2e-07 );
   petscSolveTest( 2, "../../data/meshes/3D/pyramid_tilted_4el.msh",     7.3e-06 );
   petscSolveTest( 3, "../../data/meshes/3D/regular_octahedron_8el.msh", 1.7e-06 );
   petscSolveTest( 2, "../../data/meshes/3D/cube_24el.msh",              1.4e-05 );

   return EXIT_SUCCESS;
}
