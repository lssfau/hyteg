/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/adaptiverefinement/mesh.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
// #include "hyteg/memory/MemoryAllocation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::uint_t;

#define R_min 1.0
#define R_max 2.0

// functions corresponding to a pde of type -∇⋅(k∇u) = f
struct PDE_data
{
   std::function< real_t( const hyteg::Point3D& ) > u_anal; // analytic solution u
   std::function< real_t( const hyteg::Point3D& ) > f;      // rhs f
   std::function< real_t( const hyteg::Point3D& ) > k;      // diffusion coefficient k
   std::function< real_t( const hyteg::Point3D& ) > u_D;    // dirichlet data
};

PDE_data functions( uint_t dim, uint_t shape, real_t alpha, real_t beta )
{
   PDE_data pde;

   pde.u_anal = [=]( const hyteg::Point3D& x ) {
      auto x0 = tanh( alpha * ( x.norm() - 1.5 ) );
      return 0.5 * ( exp( 2 * beta ) * std::expint( -beta * ( x0 + 1 ) ) - std::expint( -beta * ( x0 - 1 ) ) ) * exp( -beta ) /
             alpha;
   };

   pde.k = [=]( const hyteg::Point3D& x ) {
      auto x0 = x.norm();
      auto x1 = ( dim == 2 ) ? x0 : x0 * x0;
      return exp( beta * tanh( alpha * ( x0 - 1.5 ) ) ) / x1;
   };

   pde.f = [=]( const hyteg::Point3D& ) { return 0; };

   if ( shape == 0 ) // no blending
   {
      pde.u_D = pde.u_anal;
   }
   else // blending
   {
      auto R_mid   = ( R_min + R_max ) / 2.0;
      auto u_inner = pde.u_anal( Point3D( { R_min, 0, 0 } ) );
      auto u_outer = pde.u_anal( Point3D( { R_max, 0, 0 } ) );

      pde.u_D = [=]( const hyteg::Point3D& x ) { return ( x.norm() < R_mid ) ? u_inner : u_outer; };
   }

   return pde;
}

SetupPrimitiveStorage domain( uint_t dim, uint_t shape, uint_t N1, uint_t N2, uint_t N3 )
{
   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();

   if ( dim == 3 && shape == 0 )
   {
      Point3D n( { 1, 1, 1 } );
      n /= n.norm();
      meshInfo = MeshInfo::meshCuboid( R_min * n, R_max * n, N1, N2, N3 );
   }
   else if ( dim == 3 && shape == 1 )
   {
      meshInfo = MeshInfo::meshSphericalShell( N1, N2, R_min, R_max );
   }

   else if ( dim == 2 && shape == 0 )
   {
      Point2D n( { 1, 1 } );
      n /= n.norm();
      meshInfo = MeshInfo::meshRectangle( R_min * n, R_max * n, MeshInfo::CRISS, N1, N2 );
   }
   else if ( dim == 2 && shape == 1 )
   {
      meshInfo = MeshInfo::meshAnnulus( R_min, R_max, MeshInfo::CRISS, N1, N2 );
   }
   else
   {
      WALBERLA_ABORT( "Dimension must be either 2 or 3, shape must be either 0 or 1!" );
   }

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // apply geometry map
   if ( shape == 1 )
   {
      if ( dim == 3 )
      {
         IcosahedralShellMap::setMap( setupStorage );
      }
      else // dim == 2
      {
         AnnulusMap::setMap( setupStorage );
      }
   }

   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   return setupStorage;
}

// solve problem with current refinement and return sorted list of elementwise squared errors
std::vector< std::pair< real_t, hyteg::PrimitiveID > > solve( std::shared_ptr< PrimitiveStorage > storage,
                                                              const PDE_data&                     pde,
                                                              uint_t                              l_min,
                                                              uint_t                              l_max,
                                                              uint_t                              max_iter,
                                                              real_t                              tol,
                                                              std::string                         vtkname,
                                                              uint_t                              refinement_step,
                                                              bool                                l2_error_each_iteration = true )
{
   uint_t dim = storage->hasGlobalCells() ? 3 : 2;

   // operators
   using M_t = P1BlendingMassOperator;
   using A_t = P1BlendingDivkGradOperator;
   using P_t = P1toP1LinearProlongation;
   using R_t = P1toP1LinearRestriction;

   forms::p1_div_k_grad_blending_q3 form( pde.k, pde.k );

   auto M = std::make_shared< M_t >( storage, l_min, l_max + 1 );
   auto A = std::make_shared< A_t >( storage, l_min, l_max + 1, form );
   auto P = std::make_shared< P_t >();
   auto R = std::make_shared< R_t >();

   // FE functions
   P1Function< real_t > b( "b", storage, l_min, l_max );
   P1Function< real_t > u( "u", storage, l_min, l_max + 1 );
   P1Function< real_t > r( "r", storage, l_max, l_max );
   P1Function< real_t > u_anal( "u_anal", storage, l_max, l_max + 1 );
   P1Function< real_t > err( "err", storage, l_max, l_max + 1 );
   P1Function< real_t > tmp( "tmp", storage, l_min, l_max + 1 );

   // global DoF
   tmp.interpolate( []( const hyteg::Point3D& ) { return 1.0; }, l_max, hyteg::Inner );
   auto n_dof = uint_t( tmp.dotGlobal( tmp, l_max ) );
   WALBERLA_LOG_INFO_ON_ROOT( " -> number of global DoF: " << n_dof );

   // rhs
   tmp.interpolate( pde.f, l_max );
   M->apply( tmp, b, l_max, hyteg::All );
   // analytical solution
   u_anal.interpolate( pde.u_anal, l_max + 1 );
   // initialize u
   u.setToZero( l_max );
   u.setToZero( l_max + 1 );
   u.interpolate( pde.u_D, l_max, hyteg::DirichletBoundary );
   u.interpolate( pde.u_D, l_max + 1, hyteg::DirichletBoundary );

   // solver
   tmp.interpolate( []( const hyteg::Point3D& ) { return 1.0; }, l_min, hyteg::Inner );
   // cg_iter >= n coarse grid DoFs -> "exact" solve
   auto cg_iter  = std::max( max_iter, uint_t( tmp.dotGlobal( tmp, l_min ) ) );
   auto cg       = std::make_shared< CGSolver< A_t > >( storage, l_min, l_max, cg_iter, tol / 10 );
   auto smoother = std::make_shared< GaussSeidelSmoother< A_t > >();

   GeometricMultigridSolver< A_t > gmg( storage, smoother, cg, R, P, l_min, l_max, 3, 3 );

   // computation of residual and L2 error
   err.setToZero( l_max );
   err.setToZero( l_max + 1 );
   tmp.setToZero( l_max );
   tmp.setToZero( l_max + 1 );
   auto compute_residual = [&]() -> real_t {
      A->apply( u, tmp, l_max, hyteg::Inner, Replace );
      r.assign( { 1.0, -1.0 }, { b, tmp }, l_max, hyteg::Inner );
      return std::sqrt( r.dotGlobal( r, l_max, hyteg::Inner ) );
   };
   auto compute_L2error = [&]() -> real_t {
      P->prolongate( u, l_max, hyteg::Inner );
      err.assign( { 1.0, -1.0 }, { u, u_anal }, l_max + 1, hyteg::Inner );
      M->apply( err, tmp, l_max + 1, hyteg::All, Replace );
      return std::sqrt( err.dotGlobal( tmp, l_max + 1 ) );
   };

   // solve
   WALBERLA_LOG_INFO_ON_ROOT( " -> run multigrid solver" );
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  %6s |%18s |%13s", "k", "||r_k||/||r_0||", "L2 error" ) );

   real_t norm_r0 = compute_residual();
   real_t norm_r  = 0;
   real_t l2err   = 0;
   uint_t iter    = 0;
   do
   {
      ++iter;

      gmg.solve( *A, u, b, l_max );

      norm_r = compute_residual();

      if ( l2_error_each_iteration )
      {
         l2err = compute_L2error();
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  %6d |%18.3e |%13.3e", iter, norm_r / norm_r0, l2err ) );
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  %6d |%18.3e |", iter, norm_r / norm_r0 ) );
      }

   } while ( norm_r / norm_r0 > tol && iter < max_iter );

   if ( !l2_error_each_iteration )
   {
      l2err = compute_L2error();
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( " ->  %6s  %18s |%13.3e", "", "", l2err ) );
   }

   // compute elementwise error
   std::vector< std::pair< real_t, hyteg::PrimitiveID > > err_2_elwise_loc;

   if ( dim == 3 )
   {
      for ( auto& [id, cell] : storage->getCells() )
      {
         real_t err_2_cell = vertexdof::macrocell::dot< real_t >( l_max + 1, *cell, err.getCellDataID(), err.getCellDataID(), 0 );

         // scale squared error by cell-volume
         std::array< Point3D, 3 + 1 > vertices;
         uint_t                       i = 0;
         for ( auto& vid : cell->neighborVertices() )
         {
            vertices[i] = storage->getVertex( vid )->getCoordinates();
            ++i;
         }
         err_2_cell *= adaptiveRefinement::Simplex3::volume( vertices );

         err_2_elwise_loc.push_back( { err_2_cell, id } );
      }
   }
   else // dim == 2
   {
      for ( auto& [id, face] : storage->getFaces() )
      {
         real_t err_2_face = vertexdof::macroface::dot< real_t >( l_max + 1, *face, err.getFaceDataID(), err.getFaceDataID(), 0 );

         // scale squared error by face-volume
         std::array< Point3D, 2 + 1 > vertices;
         uint_t                       i = 0;
         for ( auto& vid : face->neighborVertices() )
         {
            vertices[i] = storage->getVertex( vid )->getCoordinates();
            ++i;
         }
         err_2_face *= adaptiveRefinement::Simplex2::volume( vertices );

         err_2_elwise_loc.push_back( { err_2_face, id } );
      }
   }

   // communication
   std::vector< std::pair< real_t, hyteg::PrimitiveID > > err_2_elwise;
   std::vector< std::pair< real_t, hyteg::PrimitiveID > > err_2_elwise_other;

   walberla::mpi::SendBuffer send;
   walberla::mpi::RecvBuffer recv;

   send << err_2_elwise_loc;
   walberla::mpi::allGathervBuffer( send, recv );
   for ( int rnk = 0; rnk < walberla::mpi::MPIManager::instance()->numProcesses(); ++rnk )
   {
      recv >> err_2_elwise_other;
      err_2_elwise.insert( err_2_elwise.end(), err_2_elwise_other.begin(), err_2_elwise_other.end() );
   }

   // sort by errors
   std::sort( err_2_elwise.begin(), err_2_elwise.end() );

   // export to vtk
   if ( vtkname != "" )
   {
      P1Function< real_t > k( "k", storage, l_min, l_max );
      k.interpolate( pde.k, l_max );

      P1Function< real_t > boundary( "boundary", storage, l_min, l_max );
      boundary.interpolate( []( const hyteg::Point3D& ) { return 0.0; }, l_max, hyteg::All );
      boundary.interpolate( []( const hyteg::Point3D& ) { return 1.0; }, l_max, hyteg::DirichletBoundary );

      u_anal.interpolate( pde.u_anal, l_max );

      R->restrict( err, l_max + 1, hyteg::Inner );

      P1Function< real_t > err_2( "err^2", storage, l_min, l_max );
      err_2.multElementwise( { err, err }, l_max, hyteg::All );

      VTKOutput vtkOutput( "output", vtkname, storage );
      vtkOutput.setVTKDataFormat( vtk::DataFormat::BINARY );
      vtkOutput.add( u );
      vtkOutput.add( k );
      vtkOutput.add( err );
      vtkOutput.add( err_2 );
      vtkOutput.add( u_anal );
      vtkOutput.add( b );
      vtkOutput.add( boundary );
      vtkOutput.write( l_max, refinement_step );
   }

   return err_2_elwise;
}

void solve_for_each_refinement( const SetupPrimitiveStorage& initialStorage,
                                const PDE_data&              pde,
                                uint_t                       n_ref,
                                uint_t                       n_el_max,
                                real_t                       p_ref,
                                uint_t                       l_min,
                                uint_t                       l_max,
                                uint_t                       max_iter,
                                real_t                       tol,
                                std::string                  vtkname )
{
   // construct adaptive mesh
   adaptiveRefinement::Mesh mesh( initialStorage );

   uint_t refinement = 0;
   while ( 1 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "* solving system with " << mesh.n_elements() << " macro elements ..." );

      // solve for current refinement
      auto storage      = mesh.make_storage();
      auto local_errors = solve( storage, pde, l_min, l_max, max_iter, tol, vtkname, refinement );

      // stop loop when reaching maximum refinement
      auto N_tot = mesh.n_elements();
      if ( refinement >= n_ref || N_tot >= n_el_max )
      {
         break;
      }

      ++refinement;
      WALBERLA_LOG_INFO_ON_ROOT( "* refinement " << refinement );
      WALBERLA_LOG_INFO_ON_ROOT( " -> n_el_old = " << N_tot);

      // compute number of elements to refine
      // real_t growth_est = (mesh.dim() == 2)? 3 : 7; // only considering red refinement
      // auto N_ref_max = real_t( n_el_max - N_tot ) / growth_est;
      // auto N_ref     = uint_t( std::ceil( std::min( real_t( N_tot ) * p_ref, N_ref_max ) ) );

      auto N_ref = uint_t( std::ceil( real_t( N_tot ) * p_ref ) );

      // WALBERLA_LOG_INFO_ON_ROOT( " -> " << N_ref << " of " << N_tot << " elements are being refined ..." );

      // collect ids of elements for refinement
      std::vector< PrimitiveID > R( N_ref );
      for ( uint_t i = 0; i < N_ref; ++i )
      {
         R[i] = local_errors[N_tot - 1 - i].second;
      }

      // apply refinement
      auto ref = mesh.refineRG( R, n_el_max );

      WALBERLA_LOG_INFO_ON_ROOT( " -> unrefined: " << ref.n_U << ", red: " << ref.n_R << ", green: " << ref.n_G << " -> n_el_new = " << ref.n_el());

      if ( ref.n_R == 0 )
      {
         break;
      }
   }
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   walberla::shared_ptr< walberla::config::Config > cfg;

   if ( argc == 1 )
   {
      walberla::shared_ptr< walberla::config::Config > cfg_( new walberla::config::Config );
      cfg_->readParameterFile( "../../hyteg/data/param/adaptiveRefinement.prm" );
      cfg = cfg_;
   }
   else
   {
      cfg = walberla::config::create( argc, argv );
   }

   // read parameters
   WALBERLA_LOG_INFO_ON_ROOT( "config = " << *cfg );
   walberla::Config::BlockHandle parameters = cfg->getOneBlock( "Parameters" );

   const uint_t dim   = parameters.getParameter< uint_t >( "dim" );
   const uint_t shape = parameters.getParameter< uint_t >( "shape" );
   const uint_t N1    = parameters.getParameter< uint_t >( "n1" );
   const uint_t N2    = parameters.getParameter< uint_t >( "n2" );
   const uint_t N3    = parameters.getParameter< uint_t >( "n3" );

   const real_t alpha = parameters.getParameter< real_t >( "alpha" );
   const real_t beta  = parameters.getParameter< real_t >( "beta" );

   const uint_t n_refinements = parameters.getParameter< uint_t >( "n_refinements" );
   const uint_t n_el_max      = parameters.getParameter< uint_t >( "n_el_max" );
   const real_t p_refinement  = parameters.getParameter< real_t >( "n_el_ref_relative" );
   const uint_t l_max         = parameters.getParameter< uint_t >( "microlevel" );
   const uint_t l_min         = ( l_max < 2 ) ? l_max : 2;

   const uint_t max_iter = parameters.getParameter< uint_t >( "n_iterations" );
   const real_t tol      = parameters.getParameter< real_t >( "tolerance" );

   const std::string vtkname = parameters.getParameter< std::string >( "vtkName" );

   // solve
   auto setupStorage = domain( dim, shape, N1, N2, N3 );
   auto pde_data     = functions( dim, shape, alpha, beta );
   solve_for_each_refinement(
       setupStorage, pde_data, n_refinements, n_el_max, p_refinement, l_min, l_max, max_iter, tol, vtkname );

   return 0;
}
