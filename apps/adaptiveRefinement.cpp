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
#include <core/timing/Timer.h>

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/mesh/adaptiverefinement/mesh.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using namespace hyteg;
using walberla::real_t;
using walberla::uint_t;

#define PI 3.14159265359

// solve problem with current refinement and return sorted list of elementwise squared l2 errors
std::vector< std::pair< real_t, hyteg::PrimitiveID > >
    solve( std::shared_ptr< PrimitiveStorage > storage, uint_t iter, real_t tol, const std::string& vtk )
{
   uint_t l_min = 0;
   uint_t l_max = l_min;

   // auto coeff = [=]( const hyteg::Point3D& ) { return 1; };
   auto exact = [=]( const hyteg::Point3D& x ) { return tanh( sqrt( pow( x[0], 2 ) + pow( x[1], 2 ) + pow( x[2], 2 ) ) ); };
   auto rhs   = [=]( const hyteg::Point3D& x ) {
      auto x0 = pow( x[0], 2 );
      auto x1 = pow( x[1], 2 );
      auto x2 = pow( x[2], 2 );
      auto x3 = x0 + x1 + x2;
      auto x4 = sqrt( x3 );
      auto x5 = tanh( x4 );
      auto x6 = 1 - pow( x5, 2 );
      auto x7 = x6 / pow( x3, 3.0 / 2.0 );
      auto x8 = 2 * x5 * x6 / x3;
      return x0 * x7 + x0 * x8 + x1 * x7 + x1 * x8 + x2 * x7 + x2 * x8 - 3 * x6 / x4;
   };

   P1Function< real_t > b( "b", storage, l_min, l_max );
   P1Function< real_t > u( "u", storage, l_min, l_max );
   P1Function< real_t > u_exact( "u_exact", storage, l_min, l_max );
   P1Function< real_t > err( "err", storage, l_min, l_max );
   P1Function< real_t > tmp( "tmp", storage, l_min, l_max );

   // rhs
   P1BlendingMassOperator M( storage, l_min, l_max );
   tmp.interpolate( rhs, l_max );
   M.apply( tmp, b, l_max, hyteg::All );

   // exact solution
   u_exact.interpolate( exact, l_max );

   // initialize u
   u.interpolate( exact, l_max, hyteg::DirichletBoundary );
   u.interpolate( []( const hyteg::Point3D& ) { return 0.0; }, l_max, hyteg::Inner );

   // define operator and solver
   P1BlendingLaplaceOperator                    A( storage, l_min, l_max );
   hyteg::CGSolver< P1BlendingLaplaceOperator > solver( storage, l_min, l_min, iter, tol );

   // solve
   solver.solve( A, u, b, l_max );

   // compute total error
   err.assign( { 1.0, -1.0 }, { u, u_exact }, l_max );
   M.apply( err, tmp, l_max, hyteg::All, Replace );
   real_t l2err_tot = std::sqrt( err.dotGlobal( tmp, l_max ) );

   WALBERLA_LOG_INFO_ON_ROOT( "L2-error = " << l2err_tot );

   // compute elementwise error
   std::vector< std::pair< real_t, hyteg::PrimitiveID > > l2_err_2_elwise;
   if ( storage->hasGlobalCells() )
   {
      for ( auto& [id, cell] : storage->getCells() )
      {
         real_t l2_err_2_cell = vertexdof::macrocell::dot< real_t >( l_max, *cell, err.getCellDataID(), tmp.getCellDataID() );
         l2_err_2_elwise.push_back( { l2_err_2_cell, id } );
      }
   }
   else
   {
      for ( auto& [id, face] : storage->getFaces() )
      {
         real_t l2_err_2_face = vertexdof::macroface::dot< real_t >( l_max, *face, err.getFaceDataID(), tmp.getFaceDataID() );
         l2_err_2_elwise.push_back( { l2_err_2_face, id } );
      }
   }
   // todo communicaiton
   // sort by errors
   std::sort( l2_err_2_elwise.begin(), l2_err_2_elwise.end() );

   // export to vtk
   if ( not vtk.empty() )
   {
      VTKOutput vtkOutput( "output", vtk, storage );
      vtkOutput.setVTKDataFormat( vtk::DataFormat::BINARY );
      vtkOutput.add( u );
      vtkOutput.add( err );
      vtkOutput.add( u_exact );
      vtkOutput.add( b );
      vtkOutput.write( l_max, 0 );
   }

   return l2_err_2_elwise;
}

void solve_for_each_refinement( uint_t dim, uint_t n, real_t p, uint_t iter, real_t tol, bool vtk )
{
   MeshInfo meshInfo = MeshInfo::emptyMeshInfo();

   if ( dim == 3 )
   {
      meshInfo = MeshInfo::meshCuboid( Point3D( { -1, -1, -1 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );
   }
   else if ( dim == 2 )
   {
      meshInfo = MeshInfo::meshRectangle( Point2D( { -1, -1 } ), Point2D( { 1, 1 } ), MeshInfo::CRISS, 1, 1 );
   }
   else
   {
      WALBERLA_ABORT( "Dimension must be either 2 or 3!" );
   }

   adaptiveRefinement::Mesh mesh( meshInfo );

   std::vector< std::pair< real_t, hyteg::PrimitiveID > > local_errors;

   for ( uint_t refinement = 0; refinement < n; ++refinement )
   {
      if ( refinement > 0 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "refinement " << refinement );

         uint_t N     = local_errors.size();
         uint_t N_ref = uint_t( std::ceil( real_t( N ) * p ) );

         WALBERLA_LOG_INFO_ON_ROOT( " " << N_ref << " of " << N << " elements are being refined." );

         // collect elements to refine
         std::vector< PrimitiveID > R( N_ref );
         for ( uint_t i = 0; i < N_ref; ++i )
         {
            R[i] = local_errors[N - 1 - i].second;
         }
         // apply refinement
         mesh.refineRG( R );
      }

      loadbalancing::roundRobin( mesh.setupStorage() );
      auto storage = std::make_shared< PrimitiveStorage >( mesh.setupStorage() );

      std::string vtkname = ( vtk ) ? "adaptive_r" + std::to_string( refinement ) : "";
      local_errors        = solve( storage, iter, tol, vtkname );
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

   const uint_t dim           = parameters.getParameter< uint_t >( "dim" );
   const uint_t n_refinements = parameters.getParameter< uint_t >( "n_refinements" );
   const real_t p_refinement  = parameters.getParameter< real_t >( "proportion_of_elements_refined_per_step" );
   const uint_t iter          = parameters.getParameter< uint_t >( "n_iterations" );
   const real_t tol           = parameters.getParameter< real_t >( "tolerance" );
   const bool   vtkoutput     = parameters.getParameter< bool >( "vtkOutput" );

   solve_for_each_refinement( dim, n_refinements, p_refinement, iter, tol, vtkoutput );

   return 0;
}
