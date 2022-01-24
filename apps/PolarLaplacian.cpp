/*
 * Copyright (c) 2017-2019 Christoph Schwarzmeier, Daniel Drzisga, Dominik Thoennes, Marcus Mohr.
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
#include <core/timing/Timer.h>
#include <core/Environment.h>
#include <core/math/Constants.h>

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/PolarCoordsMap.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"

#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/forms/form_fenics_generated/p1_polar_laplacian.h"

using walberla::real_t;
using walberla::uint_t;
using walberla::uint_c;
using walberla::math::pi;

using namespace hyteg;

void solve_using_geometry_map( MeshInfo&, walberla::Config::BlockHandle& );
void solve_using_pimped_form( MeshInfo&, walberla::Config::BlockHandle& );

template < typename OperatorType >
void linear_solve( OperatorType lap, P1Function< real_t >& u,
                   P1Function< real_t >& rhs,
                   P1Function< real_t >& res,
                   std::shared_ptr< PrimitiveStorage > storage,
                   walberla::Config::BlockHandle& parameters,
                   uint_t minLevel, uint_t maxLevel );

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
    WALBERLA_ABORT( "Unknown value for 'method'! Please speak English!" );
  }

  return 0;
}


// ==========================
//  solve_using_geometry_map
// ==========================
void solve_using_geometry_map( MeshInfo& meshInfo, walberla::Config::BlockHandle& parameters ) {

  // extract steering parameters
  uint_t minLevel  = parameters.getParameter<uint_t>( "minlevel"    );
  uint_t maxLevel  = parameters.getParameter<uint_t>( "maxlevel"    );
  bool   outputVTK = parameters.getParameter<bool  >( "outputVTK"   );

  // Prepare storage and set geometry mapping
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  PolarCoordsMap::setMap( setupStorage );

  loadbalancing::roundRobin( setupStorage );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  // -----------------------
  //  Problem specification
  // -----------------------

  // This is in cartesian coordinates
  std::function<real_t(const Point3D&)> solFunc =
    [](const Point3D& x) {
    real_t m = 5.0;
    real_t rho = std::sqrt( x[0] * x[0] + x[1] * x[1] );
    real_t phi = std::atan2( x[1], x[0] ) + pi;
    return std::pow(2,m)/(std::pow(2,2*m)+1)*(std::pow(rho,m)+std::pow(rho,-m))*std::sin(m*phi);
  };
  P1Function< real_t > u_exact( "u_analytic", storage, maxLevel, maxLevel );
  u_exact.interpolate( solFunc, maxLevel );

  // Create function for numeric solution and set Dirichlet boundary conditions
  std::function<real_t(const Point3D&)> zeros = [](const Point3D&) { return 0.0; };
  P1Function< real_t > u( "u_numeric", storage, minLevel, maxLevel );
  u.interpolate( zeros, maxLevel, Inner );
  u.interpolate( solFunc, maxLevel, DirichletBoundary );

  // Specify right-hand side of problem
  P1Function< real_t > rhs( "rhs", storage, minLevel, maxLevel );
  rhs.interpolate( zeros, maxLevel, All );

  // Operator for weak-form
  P1BlendingLaplaceOperator lap( storage, minLevel, maxLevel );

  // ---------------------
  //  Solve linear system
  // ---------------------
  P1Function< real_t > res( "residual", storage, minLevel, maxLevel );
  linear_solve<P1BlendingLaplaceOperator>( lap, u, rhs, res, storage, parameters, minLevel, maxLevel );

  // -------------
  //  Postprocess
  // -------------

  // compute error and its (approximate) norms
  P1Function< real_t > error( "error", storage, maxLevel, maxLevel );
  error.assign( {1.0, -1.0}, { u_exact, u }, maxLevel, All );
  real_t npoints = static_cast<real_t>( numberOfGlobalDoFs<P1FunctionTag>( *storage, maxLevel ) );
  real_t errNorm = std::sqrt( error.dotGlobal( error, maxLevel, All ) / npoints );
  real_t maxNorm = error.getMaxMagnitude( maxLevel );

  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: L_2 norm of error = " << std::scientific << errNorm );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: max norm of error = " << std::scientific << maxNorm );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: maxLevel = " << maxLevel );
  uint_t aux = numberOfGlobalInnerDoFs<P1FunctionTag>( *storage, maxLevel );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: #DoFs = " << aux );

  // output data for visualisation
  if( outputVTK ) {
    VTKOutput vtkOutput( "../output", "polar", storage );
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
  uint_t minLevel  = parameters.getParameter<uint_t>( "minlevel"    );
  uint_t maxLevel  = parameters.getParameter<uint_t>( "maxlevel"    );
  bool   outputVTK = parameters.getParameter<bool  >( "outputVTK"   );

  // Prepare storage
  SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  loadbalancing::roundRobin( setupStorage );
  std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

  // -----------------------
  //  Problem specification
  // -----------------------

  // Describe and generate exact solution (in polar coordinates)
  std::function<real_t(const Point3D&)> solFunc =
    [](const Point3D& x) { real_t m = 5.0; return pow(2,m)/(pow(2,2*m)+1)*(pow(x[0],m)+pow(x[0],-m))*sin(m*x[1]); };
  P1Function< real_t > u_exact( "u_analytic", storage, maxLevel, maxLevel );
  u_exact.interpolate( solFunc, maxLevel );

  // Create function for numeric solution and set Dirichlet boundary conditions
  std::function<real_t(const Point3D&)> zeros = [](const Point3D&) { return 0.0; };
  P1Function< real_t > u( "u_numeric", storage, minLevel, maxLevel );
  u.interpolate( zeros, maxLevel, Inner );
  u.interpolate( solFunc, maxLevel, DirichletBoundary );

  // Specify right-hand side of problem
  P1Function< real_t > rhs( "rhs", storage, minLevel, maxLevel );
  rhs.interpolate( zeros, maxLevel, All );

  // Operator for weak-form
  myPolarLapOp lap( storage, minLevel, maxLevel );

  // ---------------------
  //  Solve linear system
  // ---------------------
  P1Function< real_t > res( "residual", storage, minLevel, maxLevel );
  linear_solve< myPolarLapOp >( lap, u, rhs, res, storage, parameters, minLevel, maxLevel );

  // -------------
  //  Postprocess
  // -------------

  // compute error and its (approximate) norms
  P1Function< real_t > error( "error", storage, maxLevel, maxLevel );
  error.assign( {1.0, -1.0}, { u_exact, u }, maxLevel, All );
  real_t npoints = static_cast<real_t>( numberOfGlobalDoFs<P1FunctionTag>( *storage, maxLevel ) );
  real_t errNorm = std::sqrt( error.dotGlobal( error, maxLevel, All ) / npoints );
  real_t maxNorm = error.getMaxMagnitude( maxLevel );

  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: L_2 norm of error = " << std::scientific << errNorm );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: max norm of error = " << std::scientific << maxNorm );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: maxLevel = " << maxLevel );
  uint_t aux = numberOfGlobalInnerDoFs<P1FunctionTag>( *storage, maxLevel );
  WALBERLA_LOG_INFO_ON_ROOT( " *** MG: #DoFs = " << aux );

  // output data for visualisation
  if( outputVTK ) {
    VTKOutput vtkOutput( "../output", "polar", storage );
    vtkOutput.add( u );
    vtkOutput.add( u_exact );
    vtkOutput.add( res );
    vtkOutput.add( error );
    vtkOutput.write( maxLevel );
  }
}


// ==============
//  linear_solve
// ==============
template < typename OperatorType >
void linear_solve( OperatorType lap, P1Function< real_t >& u, P1Function< real_t >& rhs, P1Function< real_t >& res,
                   std::shared_ptr< PrimitiveStorage > storage, walberla::Config::BlockHandle& parameters,
                   uint_t minLevel, uint_t maxLevel ) {

  // extract steering parameters
  uint_t maxCycles     = parameters.getParameter<uint_t>     ( "maxCycles"   );
  real_t mgTol         = parameters.getParameter<real_t>     ( "mgTolerance" );
  std::string solver_t = parameters.getParameter<std::string>( "solver"      );

  // Compute initial residual and its norm
  lap.apply( u, res, maxLevel, Inner );
  real_t res0 = std::sqrt( res.dotGlobal( res, maxLevel, Inner ) );
  real_t resCycle = 0.0;
  bool solverConverged = false;

  // Decide on solution approach
  std::string tag;
  uint_t cycle = 0;

  if( solver_t.compare( "MG" ) == 0 ) {

    tag = "MG";
    WALBERLA_LOG_INFO_ON_ROOT( " *** MG: initial residual = " << res0 );

    // Setup geometric MG as solver
    auto smoother = std::make_shared< GaussSeidelSmoother< OperatorType >  >();
    auto coarseGridSolver = std::make_shared< CGSolver< OperatorType > >( storage, minLevel, minLevel );
    auto restrictionOperator = std::make_shared< P1toP1LinearRestriction>();
    auto prolongationOperator = std::make_shared< P1toP1LinearProlongation >();

    auto solver = GeometricMultigridSolver< OperatorType >( storage, smoother, coarseGridSolver,
                                                            restrictionOperator, prolongationOperator,
                                                            minLevel, maxLevel, 3, 3 );

    // Run MG cycles
    for( cycle = 1; cycle <= maxCycles; ++cycle ) {
      solver.solve( lap, u, rhs, maxLevel );
      lap.apply( u, res, maxLevel, Inner );
      resCycle = std::sqrt( res.dotGlobal( res, maxLevel, Inner ) );
      if( resCycle < mgTol ) {
        solverConverged = true;
        break;
      }
      WALBERLA_LOG_INFO_ON_ROOT( " *** MG: residual = " << std::scientific << resCycle );
    }
  }

  else if( solver_t.compare( "CG" ) == 0 ) {

    tag = "CG";
    WALBERLA_LOG_INFO_ON_ROOT( " *** CG: initial residual = " << res0 );
    auto cgSolver = CGSolver< OperatorType >( storage, minLevel, maxLevel, maxCycles, mgTol );
    cgSolver.solve( lap, u, rhs, maxLevel );

  }

  else if( solver_t.compare( "GS" ) == 0 ) {

    tag = "GS";
    WALBERLA_LOG_INFO_ON_ROOT( " *** GS: initial residual = " << res0 );

    // Run Gauss-Seidel iterations
    for( cycle = 1; cycle <= maxCycles; ++cycle ) {
      lap.smooth_gs( u, rhs, maxLevel, hyteg::Inner );
      lap.apply( u, res, maxLevel, hyteg::Inner ); // inefficient, only for demonstration purposes
      res.assign( {1.0,-1.0}, {rhs, res}, maxLevel, hyteg::Inner );
      resCycle = std::sqrt( res.dotGlobal( res, maxLevel, hyteg::Inner ) );
      if( resCycle < mgTol ) {
        solverConverged = true;
        break;
      }
      WALBERLA_LOG_INFO_ON_ROOT( " *** GS: residual = " << std::scientific << resCycle );
    }
  }

  else {
    WALBERLA_ABORT( "Value '" << solver_t << "' for solver not supported!" );
  }

  // Report on results
  WALBERLA_LOG_INFO_ON_ROOT( " *** " << tag << ": converged = " << ( solverConverged ? "true" : "false" ) );
  WALBERLA_LOG_INFO_ON_ROOT( " *** " << tag << ": number of iterations = " << cycle-1 );
  WALBERLA_LOG_INFO_ON_ROOT( " *** " << tag << ": final residual = " << resCycle );

}
