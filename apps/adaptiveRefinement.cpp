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
#include <core/mpi/Broadcast.h>
#include <core/timing/Timer.h>

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/mesh/adaptiverefinement/mesh.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::uint_t;

#define PI 3.14159265359
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

// solve problem with current refinement and return sorted list of elementwise squared errors
std::vector< std::pair< real_t, hyteg::PrimitiveID > >
    solve( std::shared_ptr< PrimitiveStorage > storage, const PDE_data& pde, uint_t lvl, uint_t iter, real_t tol, int vtk )
{
   uint_t l_min = lvl;
   uint_t l_max = lvl;

   uint_t dim = storage->hasGlobalCells() ? 3 : 2;

   // operators
   using M_t = P1BlendingMassOperator;
   using A_t = P1BlendingDivkGradOperator;
   forms::p1_div_k_grad_blending_q3 form( pde.k, pde.k );

   M_t M( storage, l_min, l_max );
   A_t A( storage, l_min, l_max, form );

   // FE functions
   P1Function< real_t > b( "b", storage, l_min, l_max );
   P1Function< real_t > u( "u", storage, l_min, l_max );
   P1Function< real_t > u_exact( "u_exact", storage, l_min, l_max );
   P1Function< real_t > err( "err", storage, l_min, l_max );
   P1Function< real_t > tmp( "err*M*err", storage, l_min, l_max );

   // global DoF
   tmp.interpolate( []( const hyteg::Point3D& ) { return 1.0; }, l_max, hyteg::Inner );
   auto n_dof = uint_t( tmp.dotGlobal( tmp, l_max ) );
   WALBERLA_LOG_INFO_ON_ROOT( " -> number of global DoF: " << n_dof );

   // rhs
   tmp.interpolate( pde.f, l_max );
   M.apply( tmp, b, l_max, hyteg::All );
   // exact solution
   u_exact.interpolate( pde.u_anal, l_max );
   // initialize u
   u.interpolate( pde.u_D, l_max, hyteg::DirichletBoundary );
   u.interpolate( []( const hyteg::Point3D& ) { return 0.0; }, l_max, hyteg::Inner );

   // solve
   hyteg::CGSolver< A_t > solver( storage, l_min, l_min, iter, tol );
   solver.solve( A, u, b, l_max );

   // compute total error
   err.assign( { 1.0, -1.0 }, { u, u_exact }, l_max );
   M.apply( err, tmp, l_max, hyteg::All, Replace );
   real_t l2err = std::sqrt( err.dotGlobal( tmp, l_max ) );

   WALBERLA_LOG_INFO_ON_ROOT( "L2-error = " << l2err );

   // compute elementwise error
   std::vector< std::pair< real_t, hyteg::PrimitiveID > > err_2_elwise_loc;

   if ( dim == 3 )
   {
      for ( auto& [id, cell] : storage->getCells() )
      {
         real_t err_2_cell = vertexdof::macrocell::dot< real_t >( l_max, *cell, err.getCellDataID(), err.getCellDataID(), 0 );

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
         real_t err_2_face = vertexdof::macroface::dot< real_t >( l_max, *face, err.getFaceDataID(), err.getFaceDataID(), 0 );

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
   if ( vtk >= 0 )
   {
      // compute L2 error contribution of each DoF, i.e. tmp_i=err_i*[M*err]_i
      tmp.multElementwise( { err, tmp }, l_max );
      // sanity check
      auto check = std::sqrt( tmp.sumGlobal( l_max ) );
      if ( std::abs( check - l2err ) > 1e-10 )
      {
         WALBERLA_LOG_WARNING_ON_ROOT( "sanity check failed: l2err" << l2err << " != " << check );
      }

      P1Function< real_t > k( "k", storage, l_min, l_max );
      P1Function< real_t > boundary( "boundary", storage, l_min, l_max );
      k.interpolate( pde.k, l_max );
      boundary.interpolate( []( const hyteg::Point3D& ) { return 0.0; }, l_max, hyteg::All );
      boundary.interpolate( []( const hyteg::Point3D& ) { return 1.0; }, l_max, hyteg::DirichletBoundary );

      VTKOutput vtkOutput( "output", "adaptive_" + std::to_string( dim ) + "d", storage );
      vtkOutput.setVTKDataFormat( vtk::DataFormat::BINARY );
      vtkOutput.add( u );
      vtkOutput.add( k );
      vtkOutput.add( err );
      vtkOutput.add( tmp );
      vtkOutput.add( u_exact );
      vtkOutput.add( b );
      vtkOutput.add( boundary );
      vtkOutput.write( l_max, uint_t( vtk ) );
   }

   return err_2_elwise;
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
      std::vector< double > layers( N2, 1.0 / real_t(N2) );
      meshInfo = MeshInfo::meshSphericalShell( N1, layers );
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

   // setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   return setupStorage;
}

void solve_for_each_refinement( const SetupPrimitiveStorage& initialStorage,
                                const PDE_data&              pde,
                                uint_t                       n,
                                real_t                       p,
                                uint_t                       lvl,
                                uint_t                       iter,
                                real_t                       tol,
                                bool                         vtk )
{
   // construct adaptive mesh
   adaptiveRefinement::Mesh mesh( initialStorage );
   SetupPrimitiveStorage&   setupStorage = mesh.setupStorage();

   std::vector< std::pair< real_t, hyteg::PrimitiveID > > local_errors;

   // loop over refinement levels
   for ( uint_t refinement = 0; refinement <= n; ++refinement )
   {
      // apply refinement
      if ( refinement > 0 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "* refinement " << refinement );

         uint_t N_tot = local_errors.size();
         uint_t N_ref = uint_t( std::ceil( real_t( N_tot ) * p ) );

         WALBERLA_LOG_INFO_ON_ROOT( " -> " << N_ref << " of " << N_tot << " elements are being refined ..." );

         // collect elements to refine
         std::vector< PrimitiveID > R( N_ref );
         for ( uint_t i = 0; i < N_ref; ++i )
         {
            R[i] = local_errors[N_tot - 1 - i].second;
         }

         // apply refinement and update setupStorage
         setupStorage = mesh.refineRG( R );
      }

      WALBERLA_LOG_INFO_ON_ROOT( "* solving system with " << mesh.n_elements() << " macro elements ..." );

      // std::stringstream ss;
      // setupStorage.toStream( ss, true );
      // WALBERLA_LOG_INFO_ON_ROOT( ss.str() );

      // construct PrimitiveStorage from setupStorage corresponding to current refinement
      auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

      int vtkname  = ( vtk ) ? int( refinement ) : -1;
      local_errors = solve( storage, pde, lvl, iter, tol, vtkname );
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
   const real_t p_refinement  = parameters.getParameter< real_t >( "proportion_of_elements_refined_per_step" );
   const uint_t lvl           = parameters.getParameter< uint_t >( "microlevel" );

   const uint_t iter = parameters.getParameter< uint_t >( "n_iterations" );
   const real_t tol  = parameters.getParameter< real_t >( "tolerance" );

   const bool vtkoutput = parameters.getParameter< bool >( "vtkOutput" );

   // solve
   auto setupStorage = domain( dim, shape, N1, N2, N3 );
   auto pde_data     = functions( dim, shape, alpha, beta );
   solve_for_each_refinement( setupStorage, pde_data, n_refinements, p_refinement, lvl, iter, tol, vtkoutput );

   return 0;
}
