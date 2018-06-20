#include <core/timing/Timer.h>
#include <core/Environment.h>

#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/p1functionspace/P1ElementwiseOperator.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/solvers/GeometricMultiGrid.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/communication/Syncing.hpp"
#include "tinyhhg_core/VTKWriter.hpp"

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;
using walberla::math::PI;

using namespace hhg;

int main(int argc, char* argv[])
{

  walberla::Environment walberlaEnv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  walberla::shared_ptr<walberla::config::Config> cfg( new walberla::config::Config);
  cfg->readParameterFile( "./polar.prm" );
  walberla::Config::BlockHandle parameters = cfg->getOneBlock("Parameters");

  size_t minLevel  = parameters.getParameter<size_t>( "minlevel"    );
  size_t maxLevel  = parameters.getParameter<size_t>( "maxlevel"    );
  size_t maxCycles = parameters.getParameter<size_t>( "maxCycles"   );
  real_t mgTol     = parameters.getParameter<real_t>( "mgTolerance" );
  real_t cgTol     = parameters.getParameter<real_t>( "cgTolerance" );
  bool   outputVTK = parameters.getParameter<bool  >( "outputVTK"   );

  // Mesh generation
  MeshInfo meshInfo = MeshInfo::emptyMeshInfo();

  if( parameters.isDefined( "meshFile" ) )
    {
      // Use data from a mesh file, if given in parameter file
      std::string fname = parameters.getParameter<std::string>( "meshFile" );
      meshInfo = MeshInfo::fromGmshFile( fname );
      WALBERLA_LOG_INFO_ON_ROOT( " *** Using GMSH File '" << fname << "'" );
    }
  else
    {
      // Generate annulus mesh in polar coordinates
      real_t rmin = 1.0;
      real_t rmax = 2.0;

      Point2D cornerLL( { rmin, 0.0 } );
      Point2D cornerUR( { rmax, 2.0*PI } );
      meshInfo = MeshInfo::meshRectangle( cornerLL, cornerUR, MeshInfo::CROSS, 1, 6 );
      WALBERLA_LOG_INFO_ON_ROOT( " *** Using Inline Mesher" );
    }

  SetupPrimitiveStorage setupStorage( meshInfo,
                                      uint_c ( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  hhg::loadbalancing::roundRobin( setupStorage );

  // WALBERLA_LOG_INFO_ON_ROOT( "" << setupStorage );

  std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage, timingTree);

  // Describe and generate exact solution
  std::function<real_t(const hhg::Point3D&)> solFunc =
    [](const hhg::Point3D& x) { real_t m = 5.0; return pow(2,m)/(pow(2,2*m)+1)*(pow(x[0],m)+pow(x[0],-m))*sin(m*x[1]); };
  hhg::P1Function< real_t > u_exact( "u_analytic", storage, maxLevel, maxLevel );
  u_exact.interpolate( solFunc, maxLevel );

  // Create function for numeric solution and set Dirichlet boundary conditions
  std::function<real_t(const hhg::Point3D&)> zeros = [](const hhg::Point3D&) { return 0.0; };
  hhg::P1Function< real_t > u( "u_numeric", storage, minLevel, maxLevel );
  u.interpolate( zeros, maxLevel, hhg::Inner );
  u.interpolate( solFunc, maxLevel, hhg::DirichletBoundary );

  // Specify right-hand side of problem
  hhg::P1Function< real_t > rhs( "rhs", storage, minLevel, maxLevel );
  rhs.interpolate( zeros, maxLevel, hhg::All );

  // Operator for weak-form
  hhg::P1Function< real_t > microCoordX( "microCoordX", storage, minLevel, maxLevel );
  hhg::P1Function< real_t > microCoordY( "microCoordY", storage, minLevel, maxLevel );
  std::function< real_t( const hhg::Point3D& ) > compX = []( const hhg::Point3D& pp )
    { return pp[0]; };
  std::function< real_t( const hhg::Point3D& ) > compY = []( const hhg::Point3D& pp )
    { return pp[1]; };

  for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
    {
      microCoordX.interpolate( compX, lvl );
      microCoordY.interpolate( compY, lvl );

      communication::syncFunctionBetweenPrimitives( microCoordX, lvl );
      communication::syncFunctionBetweenPrimitives( microCoordY, lvl );
    }
  hhg::P1ElementwisePolarLaplaceOperator lap( storage, {&microCoordX,&microCoordY}, minLevel, maxLevel );

  // Setup geometric MG as solver
  typedef hhg::CGSolver<hhg::P1Function<real_t>, hhg::P1ElementwisePolarLaplaceOperator> CoarseSolver;
  std::shared_ptr<CoarseSolver> coarseSolver = std::make_shared<CoarseSolver>( storage, minLevel, minLevel );

  typedef P1toP1LinearRestriction RestrictionOperator;
  RestrictionOperator restrictionOperator;

  typedef P1toP1LinearProlongation ProlongationOperator;
  ProlongationOperator prolongationOperator;

  typedef GMultigridSolver<hhg::P1Function<real_t>, hhg::P1ElementwisePolarLaplaceOperator, CoarseSolver, RestrictionOperator, ProlongationOperator > GMGSolver;
  GMGSolver solver( storage, coarseSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel );

  // Prep residual
  hhg::P1Function< real_t > res( "residual", storage, minLevel, maxLevel );

  // ================
  //  Solution phase
  // ================

  // compute initial residual and its norm
  lap.apply( u, res, maxLevel, hhg::Inner );
  real_t res0 = std::sqrt( res.dot( res, maxLevel, hhg::Inner ) );
  real_t resCycle = 0.0;
  bool mgConverged = false;

  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: initial residual = " << res0 );

  uint_t cycle = 0;

  for( cycle = 1; cycle <= maxCycles; ++cycle ) {
    solver.solve( lap, u, rhs, res, maxLevel, cgTol, 100, hhg::Inner, GMGSolver::CycleType::VCYCLE, true );
    lap.apply( u, res, maxLevel, hhg::Inner );
    // res.assign( {1.0,-1.0}, {&rhs, &res}, maxLevel, hhg::Inner );
    resCycle = std::sqrt( res.dot( res, maxLevel, hhg::Inner ) );
    if( resCycle < mgTol ) {
      mgConverged = true;
      break;
    }
    WALBERLA_LOG_INFO_ON_ROOT( " *** MG: residual = " << std::scientific << resCycle );
  }

  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: converged = " << (mgConverged == false ? "false" : "true") );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: number of cyles = " << cycle );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: final residual = " << resCycle );


  // =============
  //  Postprocess
  // =============

  // compute error norm
  hhg::P1Function< real_t > error( "error", storage, maxLevel, maxLevel );
  // hhg::P1Function< real_t > tmp( "error", storage, maxLevel, maxLevel );
  // hhg::P1MassOperator mass( storage, maxLevel, maxLevel );

  // das sollte simpler gehen (query num_dofs?)
  hhg::P1Function< real_t > npoints_helper( "npoints_helper", storage, maxLevel, maxLevel );
  std::function<real_t(const hhg::Point3D&)> ones = [](const hhg::Point3D&) { return 1.0; };
  npoints_helper.interpolate( ones, maxLevel );
  real_t npoints = npoints_helper.dot( npoints_helper, maxLevel, hhg::All );

  error.assign( {1.0, -1.0}, { &u_exact, &u }, maxLevel, hhg::All );
  // mass.apply( error, tmp, maxLevel, hhg::All );
  // real_t errNorm = std::sqrt( error.dot( tmp, maxLevel, hhg::All ) );
  real_t errNorm = std::sqrt( error.dot( error, maxLevel, hhg::All ) / npoints );
  real_t maxNorm = error.getMaxMagnitude( maxLevel );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: L_2 norm of error = " << std::scientific << errNorm );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: max norm of error = " << std::scientific << maxNorm );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: maxLevel = " << maxLevel );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: #DoFs = " << (uint_t)npoints_helper.dot( npoints_helper, maxLevel, hhg::Inner ) );

  // output data for visualisation
  if( outputVTK ) {
    hhg::VTKOutput vtkOutput(  "../output", "polar" );
    vtkOutput.add( &u );
    vtkOutput.add( &u_exact );
    vtkOutput.add( &res );
    vtkOutput.add( &error );
    vtkOutput.write( maxLevel );
  }

  return 0;
}
