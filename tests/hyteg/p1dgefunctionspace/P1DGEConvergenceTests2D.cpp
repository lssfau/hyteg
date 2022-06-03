/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "core/DataTypes.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1DGEP0StokesFunction.hpp"
#include "hyteg/composites/P1DGEP0StokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1dgefunctionspace/P1DGEOperators.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

/// Manufactured solution for u = (sin(pi x) * sin(pi y) * sin(pi (x+y)), (sin(2 pi x) * sin(2 pi y) * sin(2 pi (x+y))
/// which has homogeneous dirichlet boundary conditions on a rectangular unit triangle.
real_t testHomogeneousDirichlet( const std::string& meshFile, const uint_t& level, bool writeVTK = false )
{
   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   P1DGEFunction< idx_t > numerator( "numerator", storage, level, level );
   P1DGELaplaceOperator   L( storage, level, level );
   P1DGEMassOperator      M( storage, level, level );

   numerator.enumerate( level );

   // solution as a lambda function
   std::function< real_t( const Point3D& p ) > u_x_expr = []( const Point3D& p ) -> real_t {
      const real_t x = p[0];
      const real_t y = p[1];
      return std::sin( M_PI * x ) * std::sin( M_PI * y ) * std::sin( M_PI * ( x + y ) );
   };

   std::function< real_t( const Point3D& p ) > u_y_expr = []( const Point3D& p ) -> real_t {
      const real_t x  = p[0];
      const real_t y  = p[1];
      const real_t x0 = 2 * x;
      const real_t x1 = 2 * y;
      return std::sin( M_PI * x0 ) * std::sin( M_PI * x1 ) * std::sin( M_PI * ( x0 + x1 ) );
   };

   // rhs as a lambda function
   std::function< real_t( const Point3D& p ) > f_x_expr = []( const Point3D& p ) -> real_t {
      const real_t x  = p[0];
      const real_t y  = p[1];
      const real_t x0 = M_PI * x;
      const real_t x1 = std::sin( x0 );
      const real_t x2 = M_PI * ( x + y );
      const real_t x3 = std::pow( M_PI, 2 );
      const real_t x4 = M_PI * y;
      const real_t x5 = x3 * std::sin( x4 );
      const real_t x6 = 2 * std::cos( x2 );
      return -( x1 * x3 * x6 * std::cos( x4 ) - 4 * x1 * x5 * std::sin( x2 ) + x5 * x6 * std::cos( x0 ) );
   };
   std::function< real_t( const Point3D& p ) > f_y_expr = []( const Point3D& p ) -> real_t {
      const real_t x  = p[0];
      const real_t y  = p[1];
      const real_t x0 = 2 * x;
      const real_t x1 = 2 * y;
      const real_t x2 = M_PI * ( x0 + x1 );
      const real_t x3 = M_PI * x0;
      const real_t x4 = std::sin( x3 );
      const real_t x5 = std::pow( M_PI, 2 );
      const real_t x6 = M_PI * x1;
      const real_t x7 = x5 * std::sin( x6 );
      const real_t x8 = 8 * std::cos( x2 );
      return -( x4 * x5 * x8 * std::cos( x6 ) - 16 * x4 * x7 * std::sin( x2 ) + x7 * x8 * std::cos( x3 ) );
   };

   P1DGEFunction< real_t > u( "u", storage, level, level );
   P1DGEFunction< real_t > f( "f", storage, level, level );
   P1DGEFunction< real_t > rhs( "rhs", storage, level, level );
   P1DGEFunction< real_t > sol( "sol", storage, level, level );
   P1DGEFunction< real_t > err( "err", storage, level, level );
   P1DGEFunction< real_t > Merr( "Merr", storage, level, level );

   sol.interpolate( { u_x_expr, u_y_expr }, level, All );
   f.interpolate( { f_x_expr, f_y_expr }, level, All );

   M.apply( f, rhs, level, All, Replace );

   rhs.interpolate( 0, level, DirichletBoundary );
   u.getConformingPart()->interpolate( 0, level, DirichletBoundary );

   PETScCGSolver< P1DGELaplaceOperator > solver( storage, level, numerator );
   solver.solve( L, u, rhs, level );

   err.assign( { 1.0, -1.0 }, { u, sol }, level );

   // calculate the error in the L2 norm
   M.apply( err, Merr, level, All, Replace );
   auto discrL2 = sqrt( err.dotGlobal( Merr, level, Inner ) );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../", "P1DGEPoisson2DHomogeneousDirichletConvergenceTest", storage );
      vtk.add( u );
      vtk.add( sol );
      vtk.add( err );
      vtk.add( f );
      vtk.add( *u.getConformingPart() );
      vtk.add( *u.getDiscontinuousPart() );
      vtk.add( *numerator.getConformingPart() );
      vtk.add( *numerator.getDiscontinuousPart() );
      vtk.write( level );
   }

   return discrL2;
}

real_t testStokesHomogeneousDirichlet( const std::string& meshFile, const uint_t& level, bool writeVTK = false )
{
   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   P1DGEMassOperator M( storage, level, level );

   auto           mass_form = std::make_shared< dg::DGMassFormP0P0 >();
   dg::DGOperator M_pressure( storage, level, level, mass_form );

   P1DGEP0StokesFunction< idx_t > numerator( "numerator", storage, level, level );

   P1DGEP0StokesOperator K( storage, level, level );

   numerator.enumerate( level );

   // solution as a lambda function
   std::function< real_t( const Point3D& p ) > u_x_expr = []( const Point3D& p ) -> real_t {
      const real_t x = p[0];
      const real_t y = p[1];
      return std::sin( M_PI * x ) * std::sin( M_PI * y ) * std::sin( M_PI * ( x + y ) );
   };

   std::function< real_t( const Point3D& p ) > u_y_expr = []( const Point3D& p ) -> real_t {
      const real_t x = p[0];
      const real_t y = p[1];
      return std::sin( M_PI * x ) * std::sin( M_PI * y ) * std::sin( M_PI * ( x + y ) );
   };

   std::function< real_t( const Point3D& p ) > p_expr = []( const Point3D& p ) -> real_t {
      const real_t x = p[0];
      const real_t y = p[1];
      WALBERLA_UNUSED( x );
      WALBERLA_UNUSED( y );
      return 2 * x - y + 4;
   };

   // rhs as a lambda function
   std::function< real_t( const Point3D& p ) > f_x_expr = []( const Point3D& p ) -> real_t {
      const real_t x  = p[0];
      const real_t y  = p[1];
      const real_t x0 = M_PI * x;
      const real_t x1 = std::sin( x0 );
      const real_t x2 = M_PI * ( x + y );
      const real_t x3 = std::pow( M_PI, 2 );
      const real_t x4 = M_PI * y;
      const real_t x5 = x3 * std::sin( x4 );
      const real_t x6 = 2 * std::cos( x2 );
      return -x1 * x3 * x6 * std::cos( x4 ) + 4 * x1 * x5 * std::sin( x2 ) - x5 * x6 * std::cos( x0 ) + 2;
   };
   std::function< real_t( const Point3D& p ) > f_y_expr = []( const Point3D& p ) -> real_t {
      const real_t x  = p[0];
      const real_t y  = p[1];
      const real_t x0 = std::pow( M_PI, 2 );
      const real_t x1 = M_PI * x;
      const real_t x2 = std::sin( x1 );
      const real_t x3 = M_PI * y;
      const real_t x4 = std::sin( x3 );
      const real_t x5 = M_PI * ( x + y );
      const real_t x6 = 2 * std::cos( x5 );
      return 4 * x0 * x2 * x4 * std::sin( x5 ) - x0 * x2 * x6 * std::cos( x3 ) - x0 * x4 * x6 * std::cos( x1 ) - 1;
   };
   std::function< real_t( const Point3D& p ) > g_expr = []( const Point3D& p ) -> real_t {
      const real_t x  = p[0];
      const real_t y  = p[1];
      const real_t x0 = M_PI * y;
      const real_t x1 = M_PI * x;
      const real_t x2 = std::sin( x1 );
      const real_t x3 = M_PI * ( x + y );
      const real_t x4 = M_PI * std::sin( x3 );
      const real_t x5 = std::sin( x0 );
      return -x2 * x4 * std::cos( x0 ) - 2 * M_PI * x2 * x5 * std::cos( x3 ) - x4 * x5 * std::cos( x1 );
   };

   P1DGEP0StokesFunction< real_t > u( "u", storage, level, level );
   P1DGEP0StokesFunction< real_t > f( "f", storage, level, level );
   P1DGEP0StokesFunction< real_t > rhs( "rhs", storage, level, level );
   P1DGEP0StokesFunction< real_t > sol( "sol", storage, level, level );
   P1DGEP0StokesFunction< real_t > err( "err", storage, level, level );
   P1DGEP0StokesFunction< real_t > Merr( "Merr", storage, level, level );

   auto copyBdry = []( P1DGEP0StokesFunction< real_t > fun ) {
      fun.p().setBoundaryCondition( fun.uvw().getBoundaryCondition() );
   };

   copyBdry( u );
   copyBdry( f );
   copyBdry( rhs );
   copyBdry( sol );
   copyBdry( err );
   copyBdry( Merr );

   sol.uvw().interpolate( { u_x_expr, u_y_expr }, level, All );
   sol.p().interpolate( p_expr, level, All );
   f.uvw().interpolate( { f_x_expr, f_y_expr }, level, All );
   f.p().interpolate( g_expr, level, All );

   M.apply( f.uvw(), rhs.uvw(), level, All, Replace );
   M_pressure.apply( *f.p().getDGFunction(), *rhs.p().getDGFunction(), level, All, Replace );

   // TODO: replace by minres
   PETScMinResSolver< P1DGEP0StokesOperator > solver( storage, level, numerator );
   solver.solve( K, u, rhs, level );

   err.assign( { 1.0, -1.0 }, { u, sol }, level );

   // calculate the error in the L2 norm
   M.apply( err.uvw(), Merr.uvw(), level, All, Replace );
   auto discrL2_velocity = sqrt( err.uvw().dotGlobal( Merr.uvw(), level, Inner ) );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../", "P1DGEStokes2DHomogeneousDirichletConvergenceTest", storage );
      vtk.add( u );
      vtk.add( sol );
      vtk.add( err );
      vtk.add( f );
      vtk.add( *u.uvw().getConformingPart() );
      vtk.add( *u.uvw().getDiscontinuousPart() );
      vtk.add( *numerator.uvw().getConformingPart() );
      vtk.add( *numerator.uvw().getDiscontinuousPart() );
      vtk.write( level );
   }

   return discrL2_velocity;
}

real_t testLaplaceDirichlet( const std::string& meshFile, const uint_t& level, bool writeVTK = false )
{
   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   P1DGEFunction< idx_t > numerator( "numerator", storage, level, level );
   P1DGELaplaceOperator   L( storage, level, level );
   P1DGEMassOperator      M( storage, level, level );

   P1DGELaplaceBoundaryOperator Lboundary( storage, level, level );

   numerator.enumerate( level );

   // solution as a lambda function
   std::function< real_t( const Point3D& p ) > u_x_expr = []( const Point3D& p ) -> real_t {
      const real_t x = p[0];
      const real_t y = p[1];
      return 2 * x + y;
   };

   std::function< real_t( const Point3D& p ) > u_y_expr = []( const Point3D& p ) -> real_t {
      const real_t x = p[0];
      const real_t y = p[1];
      return -y + x;
   };

   // rhs as a lambda function
   std::function< real_t( const Point3D& p ) > f_x_expr = []( const Point3D& ) -> real_t { return 0; };
   std::function< real_t( const Point3D& p ) > f_y_expr = []( const Point3D& ) -> real_t { return 0; };

   P1DGEFunction< real_t > u( "u", storage, level, level );
   P1DGEFunction< real_t > f( "f", storage, level, level );
   P1DGEFunction< real_t > rhs( "rhs", storage, level, level );
   P1DGEFunction< real_t > sol( "sol", storage, level, level );
   P1DGEFunction< real_t > solOnBoundary( "sol", storage, level, level );
   P1DGEFunction< real_t > f2( "sol", storage, level, level );
   P1DGEFunction< real_t > err( "err", storage, level, level );
   P1DGEFunction< real_t > Merr( "Merr", storage, level, level );

   sol.interpolate( { u_x_expr, u_y_expr }, level, All );
   solOnBoundary.assign( { 1 }, { sol }, level, DirichletBoundary );
   solOnBoundary.interpolate( 0, level, Inner );
   f.interpolate( { f_x_expr, f_y_expr }, level, All );
   Lboundary.apply( sol, f2, level, All, Replace );

   M.apply( f, rhs, level, All, Replace );
   rhs.assign( { 1, 1 }, { rhs, f2 }, level, All );

   u.getConformingPart()->interpolate( { u_x_expr, u_y_expr }, level, DirichletBoundary );

   PETScCGSolver< P1DGELaplaceOperator > solver( storage, level, numerator );
   solver.solve( L, u, rhs, level );

   err.assign( { 1.0, -1.0 }, { u, sol }, level );

   // calculate the error in the L2 norm
   M.apply( err, Merr, level, All, Replace );
   auto discrL2 = sqrt( err.dotGlobal( Merr, level, Inner ) );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../", "P1DGEPoisson2DDirichletConvergenceTest", storage );
      vtk.add( u );
      vtk.add( sol );
      vtk.add( err );
      vtk.add( f );
      vtk.add( *u.getConformingPart() );
      vtk.add( *u.getDiscontinuousPart() );
      vtk.add( *numerator.getConformingPart() );
      vtk.add( *numerator.getDiscontinuousPart() );
      vtk.add( *err.getConformingPart() );
      vtk.add( *err.getDiscontinuousPart() );
      vtk.write( level );
   }

   return discrL2;
}

real_t testStokesDirichlet( const std::string& meshFile, const uint_t& level, bool writeVTK = false, bool simple = false )
{
   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   P1DGEMassOperator M( storage, level, level );

   auto           mass_form = std::make_shared< dg::DGMassFormP0P0 >();
   dg::DGOperator M_pressure( storage, level, level, mass_form );

   P1DGEP0StokesFunction< idx_t > numerator( "numerator", storage, level, level );

   P1DGEP0StokesOperator K( storage, level, level );

   P1DGELaplaceBoundaryOperator Lboundary_velocity( storage, level, level );
   P1DGEToP0DivBoundaryOperator Lboundary_pressure( storage, level, level );

   numerator.enumerate( level );

   // solution as a lambda function
   std::function< real_t( const Point3D& p ) > u_x_expr = [simple]( const Point3D& p ) -> real_t {
      const real_t x = p[0];
      const real_t y = p[1];
      if ( simple )
      {
         return y;
      }
      else
      {
         return y + 2 * std::sin( M_PI * ( x + y ) ) + 4;
      }
   };

   std::function< real_t( const Point3D& p ) > u_y_expr = [simple]( const Point3D& p ) -> real_t {
      const real_t x = p[0];
      const real_t y = p[1];
      if ( simple )
      {
         return x;
      }
      else
      {
         return -x - 2 * std::sin( M_PI * ( x + y ) ) + 3;
      }
   };

   std::function< real_t( const Point3D& p ) > p_expr = [simple]( const Point3D& p ) -> real_t {
      const real_t x = p[0];
      const real_t y = p[1];
      if ( simple )
      {
         return 0;
      }
      else
      {
         return 2 * x - y + 1;
      }
   };

   // rhs as a lambda function
   std::function< real_t( const Point3D& p ) > f_x_expr = [simple]( const Point3D& p ) -> real_t {
      const real_t x = p[0];
      const real_t y = p[1];
      if ( simple )
      {
         return 0;
      }
      else
      {
         return 4 * std::pow( M_PI, 2 ) * std::sin( M_PI * ( x + y ) ) + 2;
      }
   };
   std::function< real_t( const Point3D& p ) > f_y_expr = [simple]( const Point3D& p ) -> real_t {
      const real_t x = p[0];
      const real_t y = p[1];
      if ( simple )
      {
         return 0;
      }
      else
      {
         return -4 * std::pow( M_PI, 2 ) * std::sin( M_PI * ( x + y ) ) - 1;
      }
   };
   std::function< real_t( const Point3D& p ) > g_expr = []( const Point3D& ) -> real_t { return 0.; };

   P1DGEP0StokesFunction< real_t > u( "u", storage, level, level );
   P1DGEP0StokesFunction< real_t > f( "f", storage, level, level );
   P1DGEP0StokesFunction< real_t > rhs( "rhs", storage, level, level );
   P1DGEP0StokesFunction< real_t > sol( "sol", storage, level, level );
   P1DGEP0StokesFunction< real_t > solOnBoundary( "solOnBoundary", storage, level, level );
   P1DGEP0StokesFunction< real_t > f2( "f2", storage, level, level );
   P1DGEP0StokesFunction< real_t > err( "err", storage, level, level );
   P1DGEP0StokesFunction< real_t > Merr( "Merr", storage, level, level );

   auto copyBdry = []( P1DGEP0StokesFunction< real_t > fun ) {
      fun.p().setBoundaryCondition( fun.uvw().getBoundaryCondition() );
   };

   copyBdry( u );
   copyBdry( f );
   copyBdry( rhs );
   copyBdry( sol );
   copyBdry( err );
   copyBdry( Merr );
   copyBdry( f2 );
   copyBdry( solOnBoundary );

   sol.uvw().interpolate( { u_x_expr, u_y_expr }, level, All );
   sol.p().interpolate( p_expr, level, All );
   f.uvw().interpolate( { f_x_expr, f_y_expr }, level, All );
   f.p().interpolate( g_expr, level, All );

   solOnBoundary.uvw().assign( { 1 }, { sol.uvw() }, level, DirichletBoundary );
   solOnBoundary.interpolate( 0, level, Inner );

   Lboundary_velocity.apply( solOnBoundary.uvw(), f2.uvw(), level, All, Replace );
   Lboundary_pressure.apply( solOnBoundary.uvw(), f2.p(), level, All, Replace );

   M.apply( f.uvw(), rhs.uvw(), level, All, Replace );
   M_pressure.apply( *f.p().getDGFunction(), *rhs.p().getDGFunction(), level, All, Replace );
   rhs.assign( { 1, 1 }, { rhs, f2 }, level, All );

   u.uvw().getConformingPart()->interpolate( { u_x_expr, u_y_expr }, level, DirichletBoundary );

   PETScMinResSolver< P1DGEP0StokesOperator > solver( storage, level, numerator );
   solver.solve( K, u, rhs, level );

   err.assign( { 1.0, -1.0 }, { u, sol }, level );

   // calculate the error in the L2 norm
   M.apply( err.uvw(), Merr.uvw(), level, All, Replace );
   auto discrL2_velocity = sqrt( err.uvw().dotGlobal( Merr.uvw(), level, Inner ) );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../", "P1DGEStokes2DHomogeneousDirichletConvergenceTest", storage );
      vtk.add( u );
      vtk.add( sol );
      vtk.add( err );
      vtk.add( Merr );
      vtk.add( f );
      vtk.add( *u.uvw().getConformingPart() );
      vtk.add( *u.uvw().getDiscontinuousPart() );
      vtk.add( *numerator.uvw().getConformingPart() );
      vtk.add( *numerator.uvw().getDiscontinuousPart() );
      vtk.write( level );
   }

   return discrL2_velocity;
}

void runLaplace()
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6s|%15s|%15s", "level", "error", "rate" ) );
   real_t lastError    = std::nan( "" );
   real_t currentError = std::nan( "" );
   real_t currentRate  = std::nan( "" );
   for ( uint_t level = 6; level <= 7; level++ )
   {
      lastError    = currentError;
      currentError = hyteg::testHomogeneousDirichlet( "../../data/meshes/tri_1el.msh", level, false );
      currentRate  = lastError / currentError;
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6d|%15.2e|%15.2e", level, currentError, currentRate ) );
   }

   const real_t expectedRate = 4.;
   WALBERLA_CHECK_LESS( 0.9 * expectedRate, currentRate, "unexpected rate!" );
   WALBERLA_CHECK_GREATER( 1.1 * expectedRate, currentRate, "unexpected rate!" );
}

void runStokes()
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6s|%15s|%15s", "level", "error", "rate" ) );
   real_t lastError    = std::nan( "" );
   real_t currentError = std::nan( "" );
   real_t currentRate  = std::nan( "" );
   for ( uint_t level = 5; level <= 6; level++ )
   {
      lastError    = currentError;
      currentError = hyteg::testStokesHomogeneousDirichlet( "../../data/meshes/tri_1el.msh", level, false );
      currentRate  = lastError / currentError;
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6d|%15.2e|%15.2e", level, currentError, currentRate ) );
   }

   const real_t expectedRate = 4.;
   WALBERLA_CHECK_LESS( 0.9 * expectedRate, currentRate, "unexpected rate!" );
   WALBERLA_CHECK_GREATER( 1.1 * expectedRate, currentRate, "unexpected rate!" );
}

void runStokesDirichlet()
{
   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6s|%15s|%15s", "level", "error", "rate" ) );
   real_t lastError    = std::nan( "" );
   real_t currentError = std::nan( "" );
   real_t currentRate  = std::nan( "" );
   for ( uint_t level = 5; level <= 6; level++ )
   {
      lastError    = currentError;
      currentError = hyteg::testStokesDirichlet( "../../data/meshes/tri_1el.msh", level, false );
      currentRate  = lastError / currentError;
      WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%6d|%15.2e|%15.2e", level, currentError, currentRate ) );
   }

   const real_t expectedRate = 4.;
   WALBERLA_CHECK_LESS( 0.9 * expectedRate, currentRate, "unexpected rate!" );
   WALBERLA_CHECK_GREATER( 1.1 * expectedRate, currentRate, "unexpected rate!" );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   WALBERLA_CHECK_LESS( std::abs( hyteg::testLaplaceDirichlet( "../../data/meshes/tri_1el.msh", 4, true ) ), 1e-12 );
   WALBERLA_CHECK_LESS( std::abs( hyteg::testStokesDirichlet( "../../data/meshes/tri_1el.msh", 4, true, true ) ), 1e-12 );

   hyteg::runLaplace();
   hyteg::runStokes();
   hyteg::runStokesDirichlet();

   return 0;
}
