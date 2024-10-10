/*
 * Copyright (c) 2024 Benjamin Mann.
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
#include <core/Format.hpp>
#include <core/config/Create.h>
#include <core/math/Constants.h>
#include <core/mpi/Broadcast.h>
#include <core/timing/Timer.h>
#include <iomanip>

#include "hyteg/adaptiverefinement/error_estimator.hpp"
#include "hyteg/adaptiverefinement/mesh.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/memory/MemoryAllocation.hpp"
#include "hyteg/numerictools/L2Space.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/FullMultigridSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::uint_t;

// using Mass = P1ConstantMassOperator;
// using Mass    = P1ElementwiseMassOperator;
using Mass    = operatorgeneration::P1ElementwiseMass;
using Laplace = P1ConstantLaplaceOperator;

SetupPrimitiveStorage domain( const uint_t dim )
{
   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();

   if ( dim == 2 )
   {
      auto filename = "data/LShape_6el.msh";
      meshInfo      = MeshInfo::fromGmshFile( filename );
   }
   else
   {
      Point3D n( { 1, 1, 1 } );
      meshInfo = MeshInfo::meshCuboid( -n, n, 1, 1, 1 );
   }

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // dirichlet boundary conditions on entire boundary
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   return setupStorage;
}

void initialize( std::shared_ptr< PrimitiveStorage >        storage,
                 std::function< real_t( const Point3D& ) >& u,
                 P1Function< real_t >&                      f,
                 P1Function< real_t >&                      u_h,
                 P1Function< real_t >&                      b,
                 uint_t                                     l_min,
                 uint_t                                     l_max )
{
   using walberla::math::pi;

   if ( storage->hasGlobalCells() ) // 3d cube with rhs f=dirac
   {
      auto _f = []( const Point3D& x ) -> real_t {
         return ( x.norm() < 1e-10 ) ? std::numeric_limits< real_t >::infinity() : 0.0;
      };
      u = []( const Point3D& x ) -> real_t { return 1.0 / ( 4.0 * pi * x.norm() ); };

      // b_i = ∫fφ_i = ∫δ_0 φ_i = φ_i(0)
      auto _b = [=]( const Point3D& x ) -> real_t { return ( _f( x ) > 0.0 ) ? 1.0 : 0.0; };
      for ( uint_t lvl = l_min; lvl <= l_max; ++lvl )
      {
         f.interpolate( _f, lvl );
         b.interpolate( _b, lvl );
      }
   }
   else // L-shape domain with rhs f=0
   {
      const real_t alpha = 2.0 / 3.0;

      u = [=]( const Point3D& x ) -> real_t {
         auto r   = x.norm();
         auto phi = std::atan2( x[1], x[0] );
         if ( phi < 0 )
            phi += 2 * pi;

         return std::pow( r, alpha ) * std::sin( alpha * phi );
      };

      for ( uint_t lvl = l_min; lvl <= l_max; ++lvl )
      {
         f.interpolate( 0, lvl );
         b.interpolate( 0, lvl );
      }
   }

   // initialize u_h
   for ( uint_t lvl = l_min; lvl <= l_max; ++lvl )
      u_h.interpolate( u, lvl, DirichletBoundary );
}

// solve problem with current refinement and return list of elementwise squared errors of local elements
adaptiveRefinement::ErrorVector solve( std::shared_ptr< hyteg::PrimitiveStorage > storage,
                                       uint_t                                     l_min,
                                       uint_t                                     l_max,
                                       uint_t                                     n_cycles,
                                       uint_t                                     max_iter,
                                       uint_t                                     n1,
                                       uint_t                                     n2,
                                       real_t                                     tol,
                                       real_t                                     cg_tol,
                                       std::string                                vtkname,
                                       bool                                       writePartitioning,
                                       uint_t                                     refinement_step,
                                       bool                                       error_indicator,
                                       bool                                       global_error_estimate,
                                       uint_t                                     error_freq )
{
   // timing
   real_t t0, t1;
   real_t t_init, t_residual, t_error, t_error_indicator, t_solve;

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "* setup system ..." );

   // operators
   t0     = walberla::timing::getWcTime();
   auto A = std::make_shared< Laplace >( storage, l_min, l_max );
   auto P = std::make_shared< P1toP1LinearProlongation<> >();
   auto R = std::make_shared< P1toP1LinearRestriction<> >();

   // FE functions
   P1Function< real_t > u( "u_h", storage, l_min, l_max );   // discrete solution
   P1Function< real_t > f( "f", storage, l_min, l_max );     // rhs of strong form
   P1Function< real_t > b( "b(f)", storage, l_min, l_max );  // rhs of weak form
   P1Function< real_t > r( "r", storage, l_min, l_max );     // residual
   P1Function< real_t > tmp( "tmp", storage, l_min, l_max ); // temporary vector

   std::function< real_t( const Point3D& ) > u_anal;

   // global DoF
   tmp.interpolate( 1.0, l_max, Inner | NeumannBoundary );
   auto n_dof = uint_t( tmp.dotGlobal( tmp, l_max ) );
   // cg-dof
   tmp.interpolate( 1.0, l_min, Inner | NeumannBoundary );
   auto n_dof_coarse = uint_t( tmp.dotGlobal( tmp, l_min ) );

   // initialize analytic solution, rhs and boundary values
   initialize( storage, u_anal, f, u, b, l_min, l_max );

   // computation of residual and L2 error
   t_residual = 0;
   t_error    = 0;

   L2Space< 5, P1Function< real_t > > L2( storage, l_max );
   std::map< PrimitiveID, real_t >    err_el;

   tmp.setToZero( l_max );

   auto compute_residual = [&]() -> real_t {
      auto my_t0 = walberla::timing::getWcTime();
      A->apply( u, tmp, l_max, Inner | NeumannBoundary, Replace );
      r.assign( { 1.0, -1.0 }, { b, tmp }, l_max, Inner | NeumannBoundary );
      auto norm_r = std::sqrt( r.dotGlobal( r, l_max, Inner | NeumannBoundary ) );
      auto my_t1  = walberla::timing::getWcTime();
      t_residual += my_t1 - my_t0;
      return norm_r;
   };
   auto compute_L2error = [&]() -> real_t {
      auto my_t0 = walberla::timing::getWcTime();
      err_el.clear();
      auto err = [&]( const Point3D& x, const PrimitiveID& id ) {
         real_t ux;
         u.evaluate( x, l_max, ux, 1e-5, id );
         return u_anal( x ) - ux;
      };
      auto norm_e = L2.norm( err, err_el );
      auto my_t1  = walberla::timing::getWcTime();
      t_error += my_t1 - my_t0;
      return norm_e;
   };

   // solver

   // smoother
   auto smoother = std::make_shared< GaussSeidelSmoother< Laplace > >();
   // auto   smoother = std::make_shared< WeightedJacobiSmoother< Laplace > >( storage, l_min, l_max, 0.66 );
   // coarse grid solver
   auto cgIter = std::max( uint_t( 50 ), 2 * n_dof_coarse );
   auto cgs    = std::make_shared< CGSolver< Laplace > >( storage, l_min, l_min, cgIter, cg_tol );

   // error indicator
   t_error_indicator = 0.0;
   std::unique_ptr< adaptiveRefinement::ErrorEstimator< P1Function< real_t > > > errorEstimator;
   uint_t j_max = global_error_estimate ? l_max - l_min - 1 : 0;
   if ( error_indicator || global_error_estimate )
   {
      errorEstimator = std::make_unique< adaptiveRefinement::ErrorEstimator< P1Function< real_t > > >( u, j_max );
   }

   // multigrid
   auto gmg = std::make_shared< GeometricMultigridSolver< Laplace > >( storage, smoother, cgs, R, P, l_min, l_max, n1, n2 );
   auto fmg = std::make_shared< FullMultigridSolver< Laplace > >( storage, gmg, P, l_min, l_max, n_cycles );
   if ( error_indicator || global_error_estimate )
   {
      auto callback = [&]( uint_t lvl ) {
         real_t my_t0 = walberla::timing::getWcTime();
         errorEstimator->fmg_callback()( lvl );
         real_t my_t1 = walberla::timing::getWcTime();
         t_error_indicator += my_t1 - my_t0;
      };
      fmg = std::make_shared< FullMultigridSolver< Laplace > >(
          storage, gmg, P, l_min, l_max, n_cycles, []( uint_t ) {}, callback );
   };

   t1     = walberla::timing::getWcTime();
   t_init = t1 - t0;

   WALBERLA_LOG_INFO_ON_ROOT( " -> number of global DoF: " << n_dof );

   WALBERLA_LOG_INFO_ON_ROOT( "" );
   WALBERLA_LOG_INFO_ON_ROOT( "* solve system ..." );
   WALBERLA_LOG_INFO_ON_ROOT( " -> run multigrid solver" );
   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( " ->  %10s |%17s |%6s%d%5s", "iteration", "||r||/||b||", "||e_", l_max, "||_L2" ) );

   // initial residual
   real_t norm_0 = compute_residual();
   real_t norm_r = real_t( 1.0 );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  %10d |%17.2e |", 0, norm_r ) );

   // solve system
   t_solve = 0;
   for ( uint_t iter = 1; iter <= max_iter; ++iter )
   {
      t0 = walberla::timing::getWcTime();
      if ( iter == 1 )
         fmg->solve( *A, u, b, l_max );
      else
         gmg->solve( *A, u, b, l_max );
      t1 = walberla::timing::getWcTime();
      t_solve += t1 - t0;
      norm_r = compute_residual() / norm_0;

      bool converged = norm_r <= tol;

      if ( error_freq > 0 && ( converged || iter % error_freq == 0 ) )
      {
         auto eL2 = compute_L2error();
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  %10d |%17.2e |%12.2e", iter, norm_r, eL2 ) );
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  %10d |%17.2e |", iter, norm_r ) );
      }

      if ( converged )
         break;
   }

   // apply error indicator
   adaptiveRefinement::ErrorVector err_2_elwise_loc;
   if ( error_indicator || global_error_estimate )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "* estimate error ..." );
      t0 = walberla::timing::getWcTime();
      errorEstimator->estimate();
      if ( error_indicator )
      {
         err_2_elwise_loc = errorEstimator->eta_T_sq();
      }
      t1 = walberla::timing::getWcTime();
      t_error_indicator += t1 - t0;
   }
   if ( !error_indicator )
   {
      // use actual local L2 error
      for ( auto& [id, err] : err_el )
      {
         err_2_elwise_loc.push_back( { err, id } );
      }
   }
   // print error estimate
   if ( global_error_estimate )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  global error estimate for lvl L: ||e_%d||_L2 ≈ η_j", l_max ) );
      real_t theta_min = 1.0, theta_max = 0.0;
      for ( uint_t j = 1; j <= j_max; ++j )
      {
         auto eta   = errorEstimator->eta( j );
         auto theta = errorEstimator->theta( j );
         auto q     = -log( theta ) / log( 2.0 );

         if ( j > 1 ) // for j=1, the estimates tend to be inaccurate
         {
            theta_min = std::min( theta, theta_min );
            theta_max = std::max( theta, theta_max );
         }

         WALBERLA_LOG_INFO_ON_ROOT(
             walberla::format( " ->  j=%d: η_j = %1.2e, θ_j ≈ %1.2f ⇒ ||e||_L2 ≈ O(h^%1.2f)", j, eta, theta, q ) );
      }
      auto rho = ( j_max < 3 ) ? 0.0 : theta_min / theta_max;
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  reliability: ϱ = %1.2f", rho ) );
      if ( rho < 0.9 )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( " ->  Low reliability! Above estimates may be inaccurate!" )
      }
   }

   WALBERLA_LOG_INFO_ON_ROOT( "* Time spent for ...  " );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "  %20s: %12.3e", "system setup", t_init ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "  %20s: %12.3e", "system solve", t_solve ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "  %20s: %12.3e", "computing residual", t_residual ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "  %20s: %12.3e", "computing error", t_error ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "  %20s: %12.3e", "error indicator", t_error_indicator ) );

   printCurrentMemoryUsage();

   // export to vtk
   if ( vtkname != "" )
   {
      // error
      P1Function< real_t > err( "err", storage, l_max, l_max );
      P1Function< real_t > err_abs( "err_abs", storage, l_max, l_max );

      auto _err = [&]( const Point3D& x ) {
         real_t ux;
         u.evaluate( x, l_max, ux, 1e-5 );
         return u_anal( x ) - ux;
      };
      auto _err_abs = [&]( const Point3D& x ) { return std::abs( _err( x ) ); };

      err.interpolate( _err, l_max );
      err_abs.interpolate( _err_abs, l_max );

      // write vtk file
      VTKOutput vtkOutput( "output", vtkname, storage );
      vtkOutput.setVTKDataFormat( vtk::DataFormat::BINARY );
      vtkOutput.add( f );
      vtkOutput.add( u );
      vtkOutput.add( b );
      vtkOutput.add( err );
      vtkOutput.add( err_abs );
      vtkOutput.write( l_max, refinement_step );
   }

   if ( writePartitioning )
   {
      std::map< std::string, std::map< PrimitiveID, real_t > > realData;
      std::map< std::string, std::map< PrimitiveID, real_t > > errData;
      errData["estL2error"] = err_el;

      writeDomainPartitioningVTK(
          *storage, "output", vtkname + "_partitioning_vertices_ts" + std::to_string( refinement_step ), VTK_VERTEX, realData );
      writeDomainPartitioningVTK(
          *storage, "output", vtkname + "_partitioning_edges_ts" + std::to_string( refinement_step ), VTK_LINE, realData );
      if ( storage->hasGlobalCells() )
      {
         writeDomainPartitioningVTK(
             *storage, "output", vtkname + "_partitioning_faces_ts" + std::to_string( refinement_step ), VTK_TRIANGLE, realData );
         writeDomainPartitioningVTK(
             *storage, "output", vtkname + "_partitioning_cells_ts" + std::to_string( refinement_step ), VTK_TETRA, errData );
      }
      else
      {
         writeDomainPartitioningVTK(
             *storage, "output", vtkname + "_partitioning_faces_ts" + std::to_string( refinement_step ), VTK_TRIANGLE, errData );
      }
   }

   return err_2_elwise_loc;
}

void solve_for_each_refinement( uint_t                                 dim,
                                uint_t                                 n_ref,
                                adaptiveRefinement::RefinementStrategy ref_strat,
                                real_t                                 ref_param,
                                real_t                                 cors_param,
                                uint_t                                 l_min,
                                uint_t                                 l_max,
                                uint_t                                 n_cycles,
                                uint_t                                 max_iter,
                                uint_t                                 n1,
                                uint_t                                 n2,
                                real_t                                 tol,
                                real_t                                 cg_tol,
                                std::string                            vtkname,
                                bool                                   writePartitioning,
                                bool                                   writeMeshfile,
                                bool                                   printMeshData,
                                bool                                   error_indicator,
                                bool                                   global_error_estimate,
                                uint_t                                 error_freq )
{
   // setup domain
   auto setupStorage = domain( dim );
   // create adaptive mesh
   adaptiveRefinement::Mesh mesh( setupStorage );
   // run refinement loop
   hyteg::adaptiveRefinement::ErrorVector local_errors;
   for ( uint_t refinement = 0; refinement <= n_ref; ++refinement )
   {
      // apply refinement
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "* apply refinement " << refinement );

      auto t0 = walberla::timing::getWcTime();
      if ( refinement == 0 )
      {
         mesh.refine_regular( n_ref );
      }
      else
      {
         mesh.refineRG( local_errors, ref_strat, ref_param, cors_param, true );
      }
      auto t1    = walberla::timing::getWcTime();
      auto t_ref = t1 - t0;
      WALBERLA_LOG_INFO_ON_ROOT( " -> n_el_new = " << mesh.n_elements() );

      // compute mesh quality
      auto v_mean       = mesh.volume() / real_t( mesh.n_elements() );
      auto v_minmax     = mesh.min_max_volume();
      auto a_minmax     = mesh.min_max_angle();
      auto a_meanminmax = mesh.mean_min_max_angle();
      auto h_mean       = ( mesh.dim() == 3 ) ? pow( 6 * v_mean, 1.0 / 3.0 ) : sqrt( 2 * v_mean );
      h_mean /= real_t( 1 << l_max );
      WALBERLA_LOG_INFO_ON_ROOT( " -> el volume (min, max, mean): " << v_minmax.first << ", " << v_minmax.second << ", "
                                                                    << v_mean );
      WALBERLA_LOG_INFO_ON_ROOT( " -> angle (min, max) over all elements: " << a_minmax.first << ", " << a_minmax.second );
      WALBERLA_LOG_INFO_ON_ROOT( " -> angle (min, max) mean over all elements: " << a_meanminmax.first << ", "
                                                                                 << a_meanminmax.second );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " -> h_mean = %3.3e (fine grid)", h_mean ) );

      // load balancing
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "* apply load balancing ..." );
      t0 = walberla::timing::getWcTime();
      mesh.loadbalancing( adaptiveRefinement::GREEDY, false, true, true );
      t1                   = walberla::timing::getWcTime();
      auto t_loadbalancing = t1 - t0;

      // create storage
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "* create PrimitiveStorage ..." );
      t0                      = walberla::timing::getWcTime();
      auto storage            = mesh.make_storage();
      t1                      = walberla::timing::getWcTime();
      auto t_primitivestorage = t1 - t0;

      // print timing info
      WALBERLA_LOG_INFO_ON_ROOT( "" );
      WALBERLA_LOG_INFO_ON_ROOT( "* Time spent for ...  " );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "  %20s: %12.3e", "AMR", t_ref ) );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "  %20s: %12.3e", "loadbalancing", t_loadbalancing ) );
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "  %20s: %12.3e", "creating storage", t_primitivestorage ) );

      // output mesh info
      if ( writeMeshfile )
      {
         mesh.exportMesh( "output/" + vtkname + "_mesh_" + std::to_string( refinement ) + ".msh" );
      }
      if ( printMeshData )
      {
         WALBERLA_ROOT_SECTION()
         {
            std::ofstream file( "output/" + vtkname + "_mesh_" + std::to_string( refinement ) + ".dat" );

            file << std::setw( 13 ) << "#||centroid||" << std::setw( 13 ) << "R" << std::setw( 13 ) << "V"
                 << "\n";

            if ( mesh.dim() == 2 )
            {
               for ( auto& el : mesh.get_elements2d() )
               {
                  file << std::setw( 13 ) << el->barycenter( mesh.get_vertices() ).norm() //
                       << std::setw( 13 ) << el->radius( mesh.get_vertices() )            //
                       << std::setw( 13 ) << el->volume( mesh.get_vertices() ) << "\n";
               }
            }
            if ( mesh.dim() == 3 )
            {
               for ( auto& el : mesh.get_elements3d() )
               {
                  file << std::setw( 13 ) << el->barycenter( mesh.get_vertices() ).norm() //
                       << std::setw( 13 ) << el->radius( mesh.get_vertices() )            //
                       << std::setw( 13 ) << el->volume( mesh.get_vertices() ) << "\n";
               }
            }
         }
      }

      // solve
      local_errors = solve( storage,
                            l_min,
                            l_max,
                            n_cycles,
                            max_iter,
                            n1,
                            n2,
                            tol,
                            cg_tol,
                            vtkname,
                            writePartitioning,
                            refinement,
                            error_indicator,
                            global_error_estimate,
                            error_freq );

      // print timing tree
      if ( refinement == n_ref )
      {
         auto& timingTree = *( storage->getTimingTree() );
         if ( walberla::mpi::MPIManager::instance()->rank() == 0 )
         {
            std::cout << "\nTiming tree final iteration:\n" << timingTree << "\n";
         }
      }
   }
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   walberla::shared_ptr< walberla::config::Config > cfg( new walberla::config::Config );

   if ( argc == 1 )
   {
      cfg->readParameterFile( "data/fullAMR.prm" );
   }
   else
   {
      cfg = walberla::config::create( argc, argv );
   }

   // read parameters
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

   const uint_t dim = parameters.getParameter< uint_t >( "dim", 0 );

   const uint_t n_refinements         = parameters.getParameter< uint_t >( "n_refinements", 0 );
   const uint_t ref_strat             = parameters.getParameter< uint_t >( "refinement_strategy", 0 );
   const real_t p_refinement          = parameters.getParameter< real_t >( "p_refinement", 0.0 );
   const real_t p_coarsen             = parameters.getParameter< real_t >( "p_coarsen", 0.0 );
   const bool   error_indicator       = parameters.getParameter< bool >( "error_indicator", false );
   bool         global_error_estimate = parameters.getParameter< bool >( "global_error_estimate", false );
   uint_t       l2error               = parameters.getParameter< uint_t >( "l2error", 0 );

   const uint_t l_min = parameters.getParameter< uint_t >( "cg_level", 0 );
   const uint_t l_max = parameters.getParameter< uint_t >( "fg_level", 3 );

   const uint_t max_iter = parameters.getParameter< uint_t >( "n_iterations", 1 );
   const uint_t n_cycles = parameters.getParameter< uint_t >( "n_cycles", 1 );
   const uint_t n1       = parameters.getParameter< uint_t >( "presmooth", 2 );
   const uint_t n2       = parameters.getParameter< uint_t >( "postsmooth", 2 );
   const real_t tol      = parameters.getParameter< real_t >( "tolerance", 1e-10 );
   const real_t cg_tol   = parameters.getParameter< real_t >( "cg_tolerance", 1e-10 );

   std::string vtkname           = parameters.getParameter< std::string >( "vtkName", "" );
   const bool  writePartitioning = parameters.getParameter< bool >( "writeDomainPartitioning", false );
   const bool  writeMeshfile     = parameters.getParameter< bool >( "writeMeshfile", false );
   const bool  printMeshData     = parameters.getParameter< bool >( "printMeshData", false );

   if ( error_indicator && l_max - l_min < 1 )
   {
      WALBERLA_LOG_WARNING_ON_ROOT(
          "Local error indicator requires at least 2 multigrid levels, i.e., microlevel - cg_level >= 1." )
      WALBERLA_LOG_WARNING_ON_ROOT( "Resetting --Parameters.error_indicator=0" );
      global_error_estimate = 0;
   }
   if ( global_error_estimate && l_max - l_min < 3 )
   {
      WALBERLA_LOG_WARNING_ON_ROOT(
          "Global error estimation requires at least 2 multigrid levels, i.e., microlevel - cg_level >= 3." )
      WALBERLA_LOG_WARNING_ON_ROOT( "Resetting --Parameters.global_error_estimate=0" );
      global_error_estimate = 0;
   }
   if ( l2error == 0 && !error_indicator )
   {
      WALBERLA_LOG_WARNING_ON_ROOT( "Running without error indicator requires computation of exact error." )
      WALBERLA_LOG_WARNING_ON_ROOT( "Resetting --Parameters.l2error=max_iter!" )
      l2error = max_iter;
   }
   if ( vtkname == "auto" )
   {
      vtkname = ( dim == 2 ) ? "L_shape_full_AMR" : "dirac_3d_full_AMR";
   }

   // print parameters
   WALBERLA_LOG_INFO_ON_ROOT( "Parameters:" )
   if ( dim == 2 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "model problem", "L-shape domain" ) );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "model problem", "dirac, 3d" ) );
   }
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "number of refinement steps", n_refinements ) );
   if ( adaptiveRefinement::RefinementStrategy( ref_strat ) == adaptiveRefinement::RefinementStrategy::WEIGHTED_MEAN )
   {
      WALBERLA_LOG_INFO_ON_ROOT( " mark all elements for refinement where ||e_T||^2 > mean_T ||e_T||^2" );
      WALBERLA_LOG_INFO_ON_ROOT( " mark all elements for coarsening where ||e_T||^2 < mean_T ||e_T||^2" );
      WALBERLA_LOG_INFO_ON_ROOT( "    with weighted mean μ(x) = (∑_i i^p x_i)/(∑_i i^p) for x_1 <= x_2 <= ..." );
      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( "    where p_refinement = %1.1f and p_coarsening = %1.1f", p_refinement, p_coarsen ) );
   }
   if ( adaptiveRefinement::RefinementStrategy( ref_strat ) == adaptiveRefinement::RefinementStrategy::PERCENTILE )
   {
      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( " %30s: %3.1f%%", "proportion of elements marked for refinement", p_refinement * 100.0 ) );
      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( " %30s: %3.1f%%", "proportion of elements marked for coarsening", p_coarsen * 100.0 ) );
   }
   if ( error_indicator )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: u_%d - u_%d", "error estimate for refinement", l_max - 1, l_max ) );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: u_%d - u", "error for refinement", l_max ) );
   }
   if ( global_error_estimate )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "compute global error estimate", global_error_estimate ) );
   }
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d / %d", "level (min/max)", l_min, l_max ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "max iterations", max_iter ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "V-cycles in FMG", n_cycles ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d/%d", "number of (gs) smoothing steps", n1, n2 ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %2.1e / %2.1e", "tolerance (CG/MG)", cg_tol, tol ) );

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "write vtk output", ( vtkname != "" ) ) );
   if ( vtkname != "" )
   {
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %s", "vtk name", vtkname.c_str() ) );
   }
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "write domain partitioning", writePartitioning ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "write  mesh to file", writeMeshfile ) );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " %30s: %d", "export mesh data", printMeshData ) );

   // solve
   solve_for_each_refinement( dim,
                              n_refinements,
                              adaptiveRefinement::RefinementStrategy( ref_strat ),
                              p_refinement,
                              p_coarsen,
                              l_min,
                              l_max,
                              n_cycles,
                              max_iter,
                              n1,
                              n2,
                              tol,
                              cg_tol,
                              vtkname,
                              writePartitioning,
                              writeMeshfile,
                              printMeshData,
                              error_indicator,
                              global_error_estimate,
                              l2error );
   return 0;
}
