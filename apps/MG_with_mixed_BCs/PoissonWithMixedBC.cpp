/*
 * Copyright (c) 2020 Marcus Mohr.
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
#include <core/Environment.h>
#include <core/math/Constants.h>
#include <core/timing/Timer.h>

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/forms/form_fenics_base/P1FenicsForm.hpp"
#include "hyteg/forms/form_fenics_generated/p1_polar_laplacian.h"
#include "hyteg/geometry/PolarCoordsMap.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

// We want to solve the following problems using Multigrid
//
// D = (0,1) x (0,1)
//
//  - u_xx - u_yy = 0    , for (x,y) in D
//
//  u(x,y) = 0           , for x = 0 or x = 1
//  u(x,y) = sin(k pi x) , for y = 1
//  du(x,y) / dn = 0     , for y = 0
//
//  with analytic solution u(x,y) = sin( k pi x ) cosh( k pi y ) / cosh( k pi )
//
//
//  -v_xx - v_yy = (k^2 + m^2) pi^2 cos( k pi x ) sin( m pi y )
//
//  v(x,y) = 0           , for y = 0 or y = 1
//  dv(x,y) / dn = 0     , for x = 0 or x = 1
//
//  with analytic solution v(x,y) = cos( k pi x ) sin( m pi y )

template < typename OperatorType >
void linear_solve( OperatorType                        lap,
                   typename OperatorType::srcType&     u,
                   typename OperatorType::srcType&     rhs,
                   typename OperatorType::srcType&     res,
                   std::shared_ptr< PrimitiveStorage > storage,
                   walberla::Config::BlockHandle&      parameters,
                   uint_t                              minLevel,
                   uint_t                              maxLevel );

template < typename FunctionType >
std::shared_ptr< RestrictionOperator< FunctionType > > getRestrictionOperator()
{
   WALBERLA_ABORT( "No default restriction operator for given function type available." );
}

template <>
std::shared_ptr< RestrictionOperator< P1Function< real_t > > > getRestrictionOperator()
{
   return std::make_shared< P1toP1LinearRestriction<> >();
}

template <>
std::shared_ptr< RestrictionOperator< P2Function< real_t > > > getRestrictionOperator()
{
   return std::make_shared< P2toP2QuadraticRestriction >();
}

template < typename FunctionType >
std::shared_ptr< ProlongationOperator< FunctionType > > getProlongationOperator()
{
   WALBERLA_ABORT( "No default prolongation operator for given function type available." );
}

template <>
std::shared_ptr< ProlongationOperator< P1Function< real_t > > > getProlongationOperator()
{
   return std::make_shared< P1toP1LinearProlongation<> >();
}

template <>
std::shared_ptr< ProlongationOperator< P2Function< real_t > > > getProlongationOperator()
{
   return std::make_shared< P2toP2QuadraticProlongation >();
}

int main( int argc, char* argv[] )
{
   // ---------------
   //  General Setup
   // ---------------

   // Setup enviroment
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // Read steering parameters
   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );
   cfg->readParameterFile( "./PoissonWithMixedBC.prm" );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

   // Generate unit square mesh
   WALBERLA_LOG_INFO_ON_ROOT( " *** Using Inline Mesher" );
   Point2D  cornerLL( {real_c( 0 ), real_c( 0 )} );
   Point2D  cornerUR( {real_c( 1 ), real_c( 1 )} );
   MeshInfo meshInfo = MeshInfo::meshRectangle( cornerLL, cornerUR, MeshInfo::CROSS, 2, 2 );

   // extract steering parameters
   uint_t minLevel  = parameters.getParameter< uint_t >( "minlevel" );
   uint_t maxLevel  = parameters.getParameter< uint_t >( "maxlevel" );
   bool   outputVTK = parameters.getParameter< bool >( "outputVTK" );
   uint_t freq1     = parameters.getParameter< uint_t >( "freq1" );
   uint_t freq2     = parameters.getParameter< uint_t >( "freq2" );
   uint_t freq3     = parameters.getParameter< uint_t >( "freq3" );

   // ============
   //  Problem #1
   // ============

   WALBERLA_LOG_INFO_ON_ROOT( " ===== PROBLEM #1 =====" );

   // select FE type
   using f1_t   = P1Function< real_t >;
   using lap1_t = P1ConstantLaplaceOperator;

   // storage & boundaries
   SetupPrimitiveStorage setupStorage1( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   real_t tol       = real_c( 1e-16 );
   auto   Neumann   = [tol]( auto p ) { return p[1] < tol; };
   auto   Dirichlet = [tol]( auto p ) { return p[1] >= tol; };
   setupStorage1.setMeshBoundaryFlagsByVertexLocation( 1, Dirichlet );
   setupStorage1.setMeshBoundaryFlagsByVertexLocation( 2, Neumann );
   std::shared_ptr< PrimitiveStorage > storage1 = std::make_shared< PrimitiveStorage >( setupStorage1 );

   // analytic solution
   std::function< real_t( const Point3D& ) > solFunc1 = [freq1]( const Point3D& p ) {
      real_t x    = p[0];
      real_t y    = p[1];
      real_t freq = real_c( freq1 );
      return std::sin( freq * pi * x ) * std::cosh( freq * pi * y ) / std::cosh( freq * pi );
   };

   f1_t exact1( "analytic1", storage1, maxLevel, maxLevel );
   exact1.interpolate( solFunc1, maxLevel );

   // rhs = 0 for pure Laplace problem
   f1_t rhs1( "rhs1", storage1, minLevel, maxLevel );
   rhs1.interpolate( 0, maxLevel, All );

   // create function to hold numeric solution and set Dirichlet BCs
   f1_t x1( "numeric1", storage1, minLevel, maxLevel );
   x1.interpolate( solFunc1, maxLevel, DirichletBoundary );
   WALBERLA_LOG_INFO_ON_ROOT( " *** Setup complete" );

   // we need a Laplace operator going together with the function space
   lap1_t lap1( storage1, minLevel, maxLevel );

   // solve the linear system of equations and compute the error
   WALBERLA_LOG_INFO_ON_ROOT( " *** Solving linear system" );
   f1_t res1( "residual1", storage1, minLevel, maxLevel );
   f1_t err1( "error1", storage1, maxLevel, maxLevel );
   linear_solve( lap1, x1, rhs1, res1, storage1, parameters, minLevel, maxLevel );
   err1.assign( {real_c( 1 ), real_c( -1 )}, {exact1, x1}, maxLevel, All );

   // output data for visualisation
   if ( outputVTK )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " *** Exporting data to ./output" );
      VTKOutput vtkOutput( "./output", "results1", storage1 );

      vtkOutput.add( exact1 );
      vtkOutput.add( x1 );
      vtkOutput.add( err1 );

      vtkOutput.write( maxLevel );
   }

   // ============
   //  Problem #2
   // ============
   WALBERLA_LOG_INFO_ON_ROOT( " ===== PROBLEM #2 =====" );

   // select FE type
   using f2_t    = P2Function< real_t >;
   using lap2_t  = P2ConstantLaplaceOperator;
   using mass2_t = P2ConstantMassOperator;

   // storage & boundaries
   SetupPrimitiveStorage setupStorage2( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto Neumann2   = [tol]( auto p ) { return ( p[0] < tol ) || ( p[0] > real_c( 1 ) - tol ); };
   auto Dirichlet2 = [tol]( auto p ) { return ( p[0] >= tol ) && ( p[0] < real_c( 1 ) - tol ); };
   setupStorage2.setMeshBoundaryFlagsByVertexLocation( 1, Dirichlet2 );
   setupStorage2.setMeshBoundaryFlagsByVertexLocation( 2, Neumann2 );
   std::shared_ptr< PrimitiveStorage > storage2 = std::make_shared< PrimitiveStorage >( setupStorage2 );

   // analytic solution
   std::function< real_t( const Point3D& ) > solFunc2 = [freq2, freq3]( const Point3D& p ) {
      real_t x = p[0];
      real_t y = p[1];
      real_t k = real_c( freq2 );
      real_t m = real_c( freq3 );
      return std::cos( k * pi * x ) * std::sin( m * pi * y );
   };

   f2_t exact2( "analytic2", storage2, maxLevel, maxLevel );
   exact2.interpolate( solFunc2, maxLevel );

   // rhs for Poisson problem
   std::function< real_t( const Point3D& ) > rhsFunc2 = [freq2, freq3]( const Point3D& p ) {
      real_t x = p[0];
      real_t y = p[1];
      real_t k = real_c( freq2 );
      real_t m = real_c( freq3 );
      return ( k * k + m * m ) * pi * pi * std::cos( k * pi * x ) * std::sin( m * pi * y );
   };
   f2_t aux2( "aux2", storage2, minLevel, maxLevel );
   f2_t rhs2( "rhs2", storage2, minLevel, maxLevel );
   aux2.interpolate( rhsFunc2, maxLevel, All );
   mass2_t mass2( storage2, maxLevel, maxLevel );
   mass2.apply( aux2, rhs2, maxLevel, All );

   // create function to hold numeric solution and set Dirichlet BCs
   f2_t x2( "numeric2", storage2, minLevel, maxLevel );
   x2.interpolate( solFunc2, maxLevel, DirichletBoundary );
   WALBERLA_LOG_INFO_ON_ROOT( " *** Setup complete" );

   // we need a Laplace operator going together with the function space
   lap2_t lap2( storage2, minLevel, maxLevel );

   // solve the linear system of equations and compute the error
   WALBERLA_LOG_INFO_ON_ROOT( " *** Solving linear system" );
   f2_t res2( "residual1", storage2, minLevel, maxLevel );
   f2_t err2( "error2", storage2, maxLevel, maxLevel );
   linear_solve( lap2, x2, rhs2, res2, storage2, parameters, minLevel, maxLevel );
   err2.assign( {real_c( 1 ), real_c( -1 )}, {exact2, x2}, maxLevel, All );

   // output data for visualisation
   if ( outputVTK )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " *** Exporting data to ./output" );
      VTKOutput vtkOutput( "./output", "results2", storage2 );

      vtkOutput.add( exact2 );
      vtkOutput.add( aux2 );
      vtkOutput.add( rhs2 );
      vtkOutput.add( x2 );
      vtkOutput.add( err2 );

      vtkOutput.write( maxLevel );
   }

   return 0;
}

// ==============
//  linear_solve
// ==============
template < typename OperatorType >
void linear_solve( OperatorType                        lap,
                   typename OperatorType::srcType&     u,
                   typename OperatorType::srcType&     rhs,
                   typename OperatorType::srcType&     res,
                   std::shared_ptr< PrimitiveStorage > storage,
                   walberla::Config::BlockHandle&      parameters,
                   uint_t                              minLevel,
                   uint_t                              maxLevel )
{
   // extract steering parameters
   uint_t      maxCycles = parameters.getParameter< uint_t >( "maxCycles" );
   real_t      mgTol     = parameters.getParameter< real_t >( "mgTolerance" );
   std::string solver_t  = parameters.getParameter< std::string >( "solver" );

   // Compute initial residual and its norm
   lap.apply( u, res, maxLevel, Inner | NeumannBoundary );
   res.assign( {real_c( 1 ), real_c( -1 )}, {rhs, res}, maxLevel, Inner | NeumannBoundary );
   real_t res0            = std::sqrt( res.dotGlobal( res, maxLevel, Inner | NeumannBoundary ) );
   real_t resCycle        = 0.0;
   bool   solverConverged = false;

   // Decide on solution approach
   std::string tag;
   uint_t      cycle = 0;

   if ( solver_t.compare( "MG" ) == 0 )
   {
      tag = "MG";
      WALBERLA_LOG_INFO_ON_ROOT( " *** MG: initial residual = " << res0 );

      // Setup geometric MG as solver
      auto smoother             = std::make_shared< GaussSeidelSmoother< OperatorType > >();
      auto coarseGridSolver     = std::make_shared< CGSolver< OperatorType > >( storage, minLevel, minLevel );
      auto restrictionOperator  = getRestrictionOperator< typename OperatorType::srcType >();
      auto prolongationOperator = getProlongationOperator< typename OperatorType::srcType >();

      auto solver = GeometricMultigridSolver< OperatorType >(
          storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel, 3, 3 );

      // Run MG cycles
      for ( cycle = 1; cycle <= maxCycles; ++cycle )
      {
         solver.solve( lap, u, rhs, maxLevel );
         lap.apply( u, res, maxLevel, Inner | NeumannBoundary );
         res.assign( {real_c( 1 ), real_c( -1 )}, {rhs, res}, maxLevel, Inner | NeumannBoundary );
         resCycle = std::sqrt( res.dotGlobal( res, maxLevel, Inner | NeumannBoundary ) );
         if ( resCycle < mgTol )
         {
            solverConverged = true;
            break;
         }
         WALBERLA_LOG_INFO_ON_ROOT( " *** MG: residual = " << std::scientific << resCycle );
      }
   }

   else
   {
      WALBERLA_ABORT( "Value '" << solver_t << "' for solver not supported!" );
   }

   // Report on results
   WALBERLA_LOG_INFO_ON_ROOT( " *** " << tag << ": converged = " << ( solverConverged ? "true" : "false" ) );
   WALBERLA_LOG_INFO_ON_ROOT( " *** " << tag << ": number of iterations = " << cycle - 1 );
   WALBERLA_LOG_INFO_ON_ROOT( " *** " << tag << ": final residual = " << resCycle );
}
