#include <core/timing/Timer.h>
#include <core/Environment.h>
#include <core/math/Constants.h>

#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/geometry/PolarCoordsMap.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"
#include "tinyhhg_core/p1functionspace/P1VariableOperator.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "tinyhhg_core/solvers/CGSolver.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "tinyhhg_core/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "tinyhhg_core/solvers/GeometricMultigridSolver.hpp"
#include "tinyhhg_core/solvers/GaussSeidelSmoother.hpp"

#include "tinyhhg_core/forms/form_fenics_base/P1FenicsForm.hpp"
#include "tinyhhg_core/forms/form_fenics_generated/p1_polar_laplacian.h"

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;
using walberla::math::pi;

using namespace hyteg;

void solve_using_geometry_map( MeshInfo&, walberla::Config::BlockHandle& );
void solve_using_pimped_form( MeshInfo&, walberla::Config::BlockHandle& );

// template class hyteg::P1VariableOperator< P1FenicsForm< p1_polar_laplacian_cell_integral_0_otherwise > >;
typedef P1VariableOperator< P1FenicsForm< p1_polar_laplacian_cell_integral_0_otherwise > > myPolarLapOp;

int main(int argc, char* argv[])
{

  enum { USE_GEOMETRY_MAP, USE_PIMPED_FORM, SCREWED } approach = SCREWED;

  // ---------------
  //  General Setup
  // ---------------

  // Setup enviroment
  walberla::Environment walberlaEnv( argc, argv );
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  // Read steering parameters
  walberla::shared_ptr<walberla::config::Config> cfg( new walberla::config::Config);
  cfg->readParameterFile( "../data/param/polar.prm" );
  walberla::Config::BlockHandle parameters = cfg->getOneBlock("Parameters");

  std::string method = parameters.getParameter<std::string>( "method" );
  WALBERLA_LOG_INFO_ON_ROOT( " *** Solution approach is '" << method << "'" );

  if( method.compare( "geometryMap" ) == 0 ) {
      approach = USE_GEOMETRY_MAP;
  }
  else if( method.compare( "pimpedForm" ) == 0 ) {
      approach = USE_PIMPED_FORM;
  }
  else {
      approach = SCREWED;
  }

  // Generate annulus mesh in polar coordinates
  real_t rmin = 1.0;
  real_t rmax = 2.0;

  Point2D cornerLL( { rmin, 0.0 } );
  Point2D cornerUR( { rmax, 2.0*pi } );

  MeshInfo meshInfo = MeshInfo::meshRectangle( cornerLL, cornerUR, MeshInfo::CROSS, 1, 6 );
  WALBERLA_LOG_INFO_ON_ROOT( " *** Using Inline Mesher" );

  switch( approach ) {
  case USE_GEOMETRY_MAP:
    solve_using_geometry_map( meshInfo, parameters );
    break;
  case USE_PIMPED_FORM:
    solve_using_pimped_form( meshInfo, parameters );
    break;
  default:
    WALBERLA_ABORT( "Unknown 'approach'! Please speak English!" );
  }

  return 0;
}


// ==========================
//  solve_using_geometry_map
// ==========================
void solve_using_geometry_map( MeshInfo& meshInfo, walberla::Config::BlockHandle& parameters ) {

  // extract steering parameters
  size_t minLevel  = parameters.getParameter<size_t>( "minlevel"    );
  size_t maxLevel  = parameters.getParameter<size_t>( "maxlevel"    );
  size_t maxCycles = parameters.getParameter<size_t>( "maxCycles"   );
  real_t mgTol     = parameters.getParameter<real_t>( "mgTolerance" );
  bool   outputVTK = parameters.getParameter<bool  >( "outputVTK"   );

  // Prepare storage and set geometry mapping
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

  for( auto it : setupStorage.getFaces() ) {
    setupStorage.setGeometryMap( it.second->getID(), std::make_shared< PolarCoordsMap >() );
  }

  for( auto it : setupStorage.getEdges() ) {
    setupStorage.setGeometryMap( it.second->getID(), std::make_shared< PolarCoordsMap >() );
  }

  for( auto it : setupStorage.getVertices() ) {
    setupStorage.setGeometryMap( it.second->getID(), std::make_shared< PolarCoordsMap >() );
  }
  hyteg::loadbalancing::roundRobin( setupStorage );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  // WALBERLA_LOG_INFO_ON_ROOT( "" << setupStorage );

  // -----------------------
  //  Problem specification
  // -----------------------

  // Describe and generate exact solution

  // This is in polar coordinates
  // std::function<real_t(const hyteg::Point3D&)> solFunc =
  //   [](const hyteg::Point3D& x) { real_t m = 5.0; return pow(2,m)/(pow(2,2*m)+1)*(pow(x[0],m)+pow(x[0],-m))*sin(m*x[1]); };
  // hyteg::P1Function< real_t > u_exact( "u_analytic", storage, maxLevel, maxLevel );
  // u_exact.interpolate( solFunc, maxLevel );

  // This is in cartesian coordinates
  std::function<real_t(const hyteg::Point3D&)> solFunc =
    [](const hyteg::Point3D& x) {
    real_t m = 5.0;
    real_t rho = std::sqrt( x[0] * x[0] + x[1] * x[1] );
    real_t phi = std::atan2( x[1], x[0] ) + pi;
    return std::pow(2,m)/(std::pow(2,2*m)+1)*(std::pow(rho,m)+std::pow(rho,-m))*std::sin(m*phi);
  };
  hyteg::P1Function< real_t > u_exact( "u_analytic", storage, maxLevel, maxLevel );
  u_exact.interpolate( solFunc, maxLevel );

  // Create function for numeric solution and set Dirichlet boundary conditions
  std::function<real_t(const hyteg::Point3D&)> zeros = [](const hyteg::Point3D&) { return 0.0; };
  hyteg::P1Function< real_t > u( "u_numeric", storage, minLevel, maxLevel );
  u.interpolate( zeros, maxLevel, hyteg::Inner );
  u.interpolate( solFunc, maxLevel, hyteg::DirichletBoundary );

  // Specify right-hand side of problem
  hyteg::P1Function< real_t > rhs( "rhs", storage, minLevel, maxLevel );
  rhs.interpolate( zeros, maxLevel, hyteg::All );

  // Operator for weak-form
  P1BlendingLaplaceOperator lap( storage, minLevel, maxLevel );

  // Setup geometric MG as solver
  auto smoother = std::make_shared< hyteg::GaussSeidelSmoother< hyteg::P1BlendingLaplaceOperator>  >();
  auto coarseGridSolver = std::make_shared< hyteg::CGSolver< hyteg::P1BlendingLaplaceOperator > >( storage, minLevel, minLevel );
  auto restrictionOperator = std::make_shared< hyteg::P1toP1LinearRestriction>();
  auto prolongationOperator = std::make_shared< hyteg::P1toP1LinearProlongation >();

  auto solver = hyteg::GeometricMultigridSolver< hyteg::P1BlendingLaplaceOperator >(
       storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3 );

  // Prep residual
  hyteg::P1Function< real_t > res( "residual", storage, minLevel, maxLevel );

  // ================
  //  Solution phase
  // ================

  // compute initial residual and its norm
  lap.apply( u, res, maxLevel, hyteg::Inner );
  real_t res0 = std::sqrt( res.dotGlobal( res, maxLevel, hyteg::Inner ) );
  real_t resCycle = 0.0;
  bool mgConverged = false;

  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: initial residual = " << res0 );

  uint_t cycle = 0;

  for( cycle = 1; cycle <= maxCycles; ++cycle ) {
    solver.solve( lap, u, rhs, maxLevel );
    lap.apply( u, res, maxLevel, hyteg::Inner );
    // res.assign( {1.0,-1.0}, {&rhs, &res}, maxLevel, hyteg::Inner );
    resCycle = std::sqrt( res.dotGlobal( res, maxLevel, hyteg::Inner ) );
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
  hyteg::P1Function< real_t > error( "error", storage, maxLevel, maxLevel );
  // hyteg::P1Function< real_t > tmp( "error", storage, maxLevel, maxLevel );
  // hyteg::P1MassOperator mass( storage, maxLevel, maxLevel );

  // das sollte simpler gehen (query num_dofs?)
  hyteg::P1Function< real_t > npoints_helper( "npoints_helper", storage, maxLevel, maxLevel );
  std::function<real_t(const hyteg::Point3D&)> ones = [](const hyteg::Point3D&) { return 1.0; };
  npoints_helper.interpolate( ones, maxLevel );
  real_t npoints = npoints_helper.dotGlobal( npoints_helper, maxLevel, hyteg::All );

  error.assign( {1.0, -1.0}, { u_exact, u }, maxLevel, hyteg::All );
  // // mass.apply( error, tmp, maxLevel, hyteg::All );
  // // real_t errNorm = std::sqrt( error.dotGlobal( tmp, maxLevel, hyteg::All ) );
  real_t errNorm = std::sqrt( error.dotGlobal( error, maxLevel, hyteg::All ) / npoints );
  real_t maxNorm = error.getMaxMagnitude( maxLevel );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: L_2 norm of error = " << std::scientific << errNorm );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: max norm of error = " << std::scientific << maxNorm );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: maxLevel = " << maxLevel );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: #DoFs = " << (uint_t)npoints_helper.dotGlobal( npoints_helper, maxLevel, hyteg::Inner ) );

  // output data for visualisation
  if( outputVTK ) {
     hyteg::VTKOutput vtkOutput("../output", "polar", storage);
    vtkOutput.add( u );
    vtkOutput.add( u_exact );
    vtkOutput.add( res );
    vtkOutput.add( error );
    vtkOutput.write( maxLevel );
  }
}


// =========================
//  solve_using_pimped_form
// =========================
void solve_using_pimped_form( MeshInfo& meshInfo, walberla::Config::BlockHandle& parameters ) {

  // extract steering parameters
  size_t minLevel  = parameters.getParameter<size_t>( "minlevel"    );
  size_t maxLevel  = parameters.getParameter<size_t>( "maxlevel"    );
  size_t maxCycles = parameters.getParameter<size_t>( "maxCycles"   );
  real_t mgTol     = parameters.getParameter<real_t>( "mgTolerance" );
  bool   outputVTK = parameters.getParameter<bool  >( "outputVTK"   );

  // Prepare storage
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  hyteg::loadbalancing::roundRobin( setupStorage );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  // -----------------------
  //  Problem specification
  // -----------------------

  // Describe and generate exact solution (in polar coordinates)
  std::function<real_t(const hyteg::Point3D&)> solFunc =
    [](const hyteg::Point3D& x) { real_t m = 5.0; return pow(2,m)/(pow(2,2*m)+1)*(pow(x[0],m)+pow(x[0],-m))*sin(m*x[1]); };
  hyteg::P1Function< real_t > u_exact( "u_analytic", storage, maxLevel, maxLevel );
  u_exact.interpolate( solFunc, maxLevel );

  // Create function for numeric solution and set Dirichlet boundary conditions
  std::function<real_t(const hyteg::Point3D&)> zeros = [](const hyteg::Point3D&) { return 0.0; };
  hyteg::P1Function< real_t > u( "u_numeric", storage, minLevel, maxLevel );
  u.interpolate( zeros, maxLevel, hyteg::Inner );
  u.interpolate( solFunc, maxLevel, hyteg::DirichletBoundary );

  // Specify right-hand side of problem
  hyteg::P1Function< real_t > rhs( "rhs", storage, minLevel, maxLevel );
  rhs.interpolate( zeros, maxLevel, hyteg::All );

  // Operator for weak-form
  myPolarLapOp lap( storage, minLevel, maxLevel );

  // Setup geometric MG as solver
  auto smoother = std::make_shared< hyteg::GaussSeidelSmoother< myPolarLapOp >  >();
  auto coarseGridSolver = std::make_shared< hyteg::CGSolver< myPolarLapOp > >( storage, minLevel, minLevel );
  auto restrictionOperator = std::make_shared< hyteg::P1toP1LinearRestriction>();
  auto prolongationOperator = std::make_shared< hyteg::P1toP1LinearProlongation >();

  auto solver = hyteg::GeometricMultigridSolver< myPolarLapOp >(
       storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3 );

  // Prep residual
  hyteg::P1Function< real_t > res( "residual", storage, minLevel, maxLevel );

  // ================
  //  Solution phase
  // ================

  // compute initial residual and its norm
  lap.apply( u, res, maxLevel, hyteg::Inner );
  real_t res0 = std::sqrt( res.dotGlobal( res, maxLevel, hyteg::Inner ) );
  real_t resCycle = 0.0;
  bool mgConverged = false;

  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: initial residual = " << res0 );

  uint_t cycle = 0;

  for( cycle = 1; cycle <= maxCycles; ++cycle ) {
    solver.solve( lap, u, rhs, maxLevel );
    lap.apply( u, res, maxLevel, hyteg::Inner );
    // res.assign( {1.0,-1.0}, {&rhs, &res}, maxLevel, hyteg::Inner );
    resCycle = std::sqrt( res.dotGlobal( res, maxLevel, hyteg::Inner ) );
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
  hyteg::P1Function< real_t > error( "error", storage, maxLevel, maxLevel );
  // hyteg::P1Function< real_t > tmp( "error", storage, maxLevel, maxLevel );
  // hyteg::P1MassOperator mass( storage, maxLevel, maxLevel );

  // das sollte simpler gehen (query num_dofs?)
  hyteg::P1Function< real_t > npoints_helper( "npoints_helper", storage, maxLevel, maxLevel );
  std::function<real_t(const hyteg::Point3D&)> ones = [](const hyteg::Point3D&) { return 1.0; };
  npoints_helper.interpolate( ones, maxLevel );
  real_t npoints = npoints_helper.dotGlobal( npoints_helper, maxLevel, hyteg::All );

  error.assign( {1.0, -1.0}, { u_exact, u }, maxLevel, hyteg::All );
  // // mass.apply( error, tmp, maxLevel, hyteg::All );
  // // real_t errNorm = std::sqrt( error.dotGlobal( tmp, maxLevel, hyteg::All ) );
  real_t errNorm = std::sqrt( error.dotGlobal( error, maxLevel, hyteg::All ) / npoints );
  real_t maxNorm = error.getMaxMagnitude( maxLevel );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: L_2 norm of error = " << std::scientific << errNorm );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: max norm of error = " << std::scientific << maxNorm );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: maxLevel = " << maxLevel );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: #DoFs = " << (uint_t)npoints_helper.dotGlobal( npoints_helper, maxLevel, hyteg::Inner ) );

  // output data for visualisation
  if( outputVTK ) {
     hyteg::VTKOutput vtkOutput("../output", "polar", storage);
    vtkOutput.add( u );
    vtkOutput.add( u_exact );
    vtkOutput.add( res );
    vtkOutput.add( error );
    vtkOutput.write( maxLevel );
  }
}

