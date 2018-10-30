
/**
 *
 * References:
 *
 * [1] Haines 2011: The Jeffery–Hamel similarity solution and its relation to flow in a diverging channel
 * [2] Fraenkel 1962: Laminar flow in symmetrical channels with slightly curved walls
 *                    I. On the Jeffery-Hamel solutions for flow between plane walls
 *
 */

#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/Visualization.hpp"

#include "tinyhhg_core/composites/P2P1TaylorHoodFunction.hpp"
#include "tinyhhg_core/composites/P2P1TaylorHoodStokesOperator.hpp"

#include "tinyhhg_core/gridtransferoperators/P2P1StokesToP2P1StokesProlongation.hpp"
#include "tinyhhg_core/gridtransferoperators/P2P1StokesToP2P1StokesRestriction.hpp"

#include "tinyhhg_core/solvers/preconditioners/StokesBlockDiagonalPreconditioner.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/MinresSolver.hpp"
#include "tinyhhg_core/solvers/UzawaSolver.hpp"

#include "tinyhhg_core/petsc/PETScManager.hpp"
#include "tinyhhg_core/petsc/PETScLUSolver.hpp"

#include "tinyhhg_core/VTKWriter.hpp"

namespace hhg {

using walberla::real_c;
using walberla::real_t;
using walberla::math::PI;

void jefferyHamelFlowTest()
{
  ////////////////
  // Parameters //
  ////////////////

  // All variable names are defined as in [1].

  // free parameters

  const double eta      = 2.0;
  const double alpha    = 15.0 * (PI / 180.0);
  const uint_t minLevel = 2;
  const uint_t maxLevel = 2;

  // derived parameters

  const double Rhat_o   = 1.0;
  const double Rhat_i   = Rhat_o / eta;
  const double epsilon  = 1e-8;

  ////////////
  // Domain //
  ////////////

  WALBERLA_LOG_INFO_ON_ROOT( "[JefferyHamelHBenchmark] Creating mesh..." );

  const auto meshInfo = MeshInfo::meshAnnulus( Rhat_i, Rhat_o, - alpha, alpha, MeshInfo::meshFlavour::DIAMOND, 4, 4 );
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  WALBERLA_LOG_INFO_ON_ROOT( "[JefferyHamelHBenchmark] Defining boundaries..." );

  // setting mesh boundary flags

  BoundaryCondition bc;
  const auto noSlipBCUID  = bc.createDirichletBC( "no-slip", 1 );
  const auto inflowBCUID  = bc.createDirichletBC( "inflow",  2 );
  const auto outflowBCUID = bc.createNeumannBC  ( "outflow", 3 );
  auto onInflowBoundary  = [ alpha, Rhat_i, epsilon ]( const Point3D & p ) { return std::abs( p.norm() - Rhat_i ) < epsilon && std::abs( math::toSpherical( p )[2] ) < alpha; };
  auto onOutflowBoundary = [ alpha, Rhat_o, epsilon ]( const Point3D & p ) { return std::abs( p.norm() - Rhat_o ) < epsilon && std::abs( math::toSpherical( p )[2] ) < alpha; };
  setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
  setupStorage.setMeshBoundaryFlagsByVertexLocation( 2, onInflowBoundary, false );
  setupStorage.setMeshBoundaryFlagsByVertexLocation( 3, onOutflowBoundary, false );
  setupStorage.setMeshBoundaryFlagsInner( 0, true );

  WALBERLA_LOG_INFO_ON_ROOT( "[JefferyHamelHBenchmark] Creating distributed domain..." );

  const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

  WALBERLA_LOG_INFO_ON_ROOT( "[JefferyHamelHBenchmark] Domain VTK output..." );

  writeDomainPartitioningVTK( storage, "../../output", "JefferyHamel_Domain" );

  ///////////////
  // Functions //
  ///////////////

  P2P1TaylorHoodFunction< real_t > x            ( "x",             storage, minLevel, maxLevel, bc );
  P2P1TaylorHoodFunction< real_t > Ax           ( "Ax",            storage, minLevel, maxLevel, bc );
  P2P1TaylorHoodFunction< real_t > exactSolution( "exactSolution", storage, minLevel, maxLevel, bc );
  P2P1TaylorHoodFunction< real_t > rhs          ( "rhs",           storage, minLevel, maxLevel, bc );
  P2P1TaylorHoodFunction< real_t > r            ( "r",             storage, minLevel, maxLevel, bc );

  P2P1TaylorHoodStokesOperator stokesOperator( storage, minLevel, maxLevel );

  //////////////
  // Boundary //
  //////////////

  auto inflowProfile = [ alpha ]( const Point3D & p, const uint_t & coordinate ) -> real_t {
    const auto pSpherical = math::toSpherical( p );
    const real_t phiRescaled = pSpherical[2] / alpha;
    const real_t uSpherical = ( real_c( 3 ) / ( real_c( 4 ) * alpha ) ) * ( real_c( 1 ) - phiRescaled * phiRescaled );
    return math::toCartesian( Point3D( { uSpherical, 0.5 * PI, pSpherical[2] } ) )[ coordinate ];
  };

  auto inflowProfileU = [ inflowProfile ]( const Point3D & p ) -> real_t { return inflowProfile( p, 0 ); };
  auto inflowProfileV = [ inflowProfile ]( const Point3D & p ) -> real_t { return inflowProfile( p, 1 ); };

  auto solution = [ alpha ]( const Point3D & p, const uint_t & coordinate ) -> real_t {
    const auto pSpherical = math::toSpherical( p );
    const real_t theta = pSpherical[2];
    const real_t y = theta / alpha;
    const real_t twoAlpha = real_c(2) * alpha;
    const real_t cosTwoAlpha = std::cos(twoAlpha);
    const real_t cosTwoAlphaY = std::cos(twoAlpha * y);
    const real_t sinTwoAlpha = std::sin(twoAlpha);
    const real_t velocitySpherical = (twoAlpha * (cosTwoAlphaY - cosTwoAlpha)) / (sinTwoAlpha - twoAlpha * cosTwoAlpha);
    const auto velocityCartesian = math::toCartesian( Point3D( { velocitySpherical, 0.5 * PI, theta } ) );
    return velocityCartesian[ coordinate ];
  };

  auto solutionU = [ solution ]( const Point3D & p ) -> real_t { return solution( p, 0 ); };
  auto solutionV = [ solution ]( const Point3D & p ) -> real_t { return solution( p, 1 ); };

  WALBERLA_LOG_DEVEL_ON_ROOT( "[JefferyHamelHBenchmark] Interpolating inflow..." );

  x.u.interpolate( solutionU, maxLevel, inflowBCUID );
  x.v.interpolate( solutionV, maxLevel, inflowBCUID );

  exactSolution.u.interpolate( solutionU, maxLevel, All );
  exactSolution.v.interpolate( solutionV, maxLevel, All );

  /////////
  // VTK //
  /////////

  WALBERLA_LOG_DEVEL_ON_ROOT( "[JefferyHamelHBenchmark] Setting up VTK..." );

  VTKOutput vtkOutput( "../../output", "JefferyHamel", storage );
  vtkOutput.add( &x.u );
  vtkOutput.add( &x.v );
  vtkOutput.add( &x.p );
  vtkOutput.add( &exactSolution.u );
  vtkOutput.add( &exactSolution.v );

  WALBERLA_LOG_DEVEL_ON_ROOT( "[JefferyHamelHBenchmark] Writing initial VTK..." );

  vtkOutput.write( maxLevel, 0 );

  ////////////
  // Solver //
  ////////////

  WALBERLA_LOG_DEVEL_ON_ROOT( "[JefferyHamelHBenchmark] Setting up solvers..." );

  typedef hhg::CGSolver< typename P2P1TaylorHoodFunction< real_t >::VelocityFunction_T, typename P2P1TaylorHoodStokesOperator::VelocityOperator_T > VelocityCGSolver_T;
  typedef hhg::GMultigridSolver< typename P2P1TaylorHoodFunction< real_t >::VelocityFunction_T,
  typename P2P1TaylorHoodStokesOperator::VelocityOperator_T,
  VelocityCGSolver_T,
  typename P2P1StokesToP2P1StokesRestriction::VelocityRestriction_T,
  typename P2P1StokesToP2P1StokesProlongation::VelocityProlongation_T > VelocityGMGSolver_T;
  typename P2P1StokesToP2P1StokesRestriction::VelocityRestriction_T   velocityRestriction;
  typename P2P1StokesToP2P1StokesProlongation::VelocityProlongation_T velocityProlongation;
  auto velocityCGSolver = std::make_shared< VelocityCGSolver_T >( storage, minLevel, maxLevel );
  VelocityGMGSolver_T velocityGMGSolver( storage, velocityCGSolver, velocityRestriction, velocityProlongation, minLevel, maxLevel, 2, 2 );

  typedef hhg::StokesBlockDiagonalPreconditioner< P2P1TaylorHoodFunction< real_t >,
  typename P2P1TaylorHoodStokesOperator::VelocityOperator_T,
  VelocityGMGSolver_T,
  hhg::P1LumpedInvMassOperator > Preconditioner_T;
  hhg::P1LumpedInvMassOperator lumpedInvMassOperator( storage, minLevel, maxLevel );
  Preconditioner_T preconditioner( stokesOperator.A, velocityGMGSolver, lumpedInvMassOperator, storage, minLevel, maxLevel, 0 );
  typedef hhg::MinResSolver< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodStokesOperator, Preconditioner_T > PreconditionedMinResSolver_T;
  auto preconditionedMinResSolver = PreconditionedMinResSolver_T( storage, minLevel, maxLevel, preconditioner );

  P2P1StokesToP2P1StokesRestriction restrictionOperator;
  P2P1StokesToP2P1StokesProlongation prolongationOperator;
  typedef hhg::UzawaSolver< P2P1TaylorHoodFunction< real_t >, P2P1TaylorHoodStokesOperator, PreconditionedMinResSolver_T, P2P1StokesToP2P1StokesRestriction, P2P1StokesToP2P1StokesProlongation, false > UzawaSolver_T;
  UzawaSolver_T uzawaSolver( storage, preconditionedMinResSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3, 2, 0.001 );

  WALBERLA_LOG_DEVEL_ON_ROOT( "[JefferyHamelHBenchmark] Solving..." );

#if 0
  typedef typename P2P1TaylorHoodFunction< real_t >::Tag StokesFunctionTag_T;

  stokesOperator.apply( x, Ax, maxLevel, hhg::Inner | hhg::NeumannBoundary );
  r.assign( {1.0, -1.0}, {&rhs, &Ax}, maxLevel, hhg::Inner | hhg::NeumannBoundary );
  real_t currentResidualL2 = std::sqrt( r.dotGlobal( r, maxLevel, hhg::Inner | hhg::NeumannBoundary ) ) / real_c(hhg::numberOfGlobalDoFs< StokesFunctionTag_T >( *storage, maxLevel ));
  real_t lastResidualL2    = currentResidualL2;
  WALBERLA_LOG_INFO_ON_ROOT( "[StokesFlowSolverComparison] Initial residual: " << currentResidualL2 );

  WALBERLA_LOG_INFO_ON_ROOT( "[JefferyHamelHBenchmark] iteration | residual (L2) | convergence rate |     time " )
  WALBERLA_LOG_INFO_ON_ROOT( "[JefferyHamelHBenchmark] ----------+---------------+------------------+--------- " )

  walberla::WcTimer timer;
  for ( uint_t mgIteration = 0; mgIteration < 2; mgIteration++ )
  {
    timer.start();
    uzawaSolver.solve( stokesOperator, x, rhs, r, maxLevel, 1e-16, 10000, Inner | NeumannBoundary );
    timer.end();

    lastResidualL2 = currentResidualL2;
    stokesOperator.apply( x, Ax, maxLevel, hhg::Inner | hhg::NeumannBoundary );
    r.assign( {1.0, -1.0}, {&rhs, &Ax}, maxLevel, hhg::Inner | hhg::NeumannBoundary );
    currentResidualL2 = std::sqrt( r.dotGlobal( r, maxLevel, hhg::Inner | hhg::NeumannBoundary ) ) / real_c(hhg::numberOfGlobalDoFs< StokesFunctionTag_T >( *storage, maxLevel ));


    WALBERLA_LOG_INFO_ON_ROOT( "[JefferyHamelHBenchmark] " << std::setw(9) << mgIteration << " | "
                                                               << std::setw(13) << std::scientific << currentResidualL2 << " | "
                                                               << std::setw(16) << std::scientific << currentResidualL2 / lastResidualL2 << " | "
                                                               << std::setw(7) << std::fixed << std::setprecision(3) << timer.last() )
  }
#endif

#if 0
  preconditionedMinResSolver.solve( stokesOperator, x, rhs, r, maxLevel, 1e-16, 10000, Inner | NeumannBoundary, true );
#endif

#if 1
  PETScManager petscManager;
  rhs.u.interpolate( solutionU, maxLevel, inflowBCUID );
  rhs.v.interpolate( solutionV, maxLevel, inflowBCUID );
  auto numerator = std::make_shared<hhg::P2P1TaylorHoodFunction< PetscInt > >("numerator", storage, minLevel, maxLevel);
  numerator->enumerate( maxLevel );
  const uint_t localDoFs = hhg::numberOfLocalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel );
  const uint_t globalDoFs = hhg::numberOfGlobalDoFs< P2P1TaylorHoodFunctionTag >( *storage, maxLevel );
  PETScLUSolver<real_t, hhg::P2P1TaylorHoodFunction, hhg::P2P1TaylorHoodStokesOperator> solver(numerator, localDoFs, globalDoFs);
  solver.solve( stokesOperator, x, rhs, r, maxLevel, 1e-16, 10000, hhg::Inner | hhg::NeumannBoundary, true );
#endif


  WALBERLA_LOG_DEVEL_ON_ROOT( "[JefferyHamelHBenchmark] Writing final VTK..." );

  vtkOutput.write( maxLevel, 1 );


}

}

int main( int argc, char** argv )
{
  walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
  walberla::MPIManager::instance()->useWorldComm();
  hhg::jefferyHamelFlowTest();
  return EXIT_SUCCESS;
}