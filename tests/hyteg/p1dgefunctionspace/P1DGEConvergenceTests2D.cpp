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

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGDiffusionForm_Example.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1dgefunctionspace/P1DGEOperators.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

/// Manufactured solution for u = (sin(pi x) * sin(pi y) * sin(pi (x+y)), (sin(2 pi x) * sin(2 pi y) * sin(2 pi (x+y))
/// which has homogeneous dirichlet boundary conditions on a rectangular unit triangle.
real_t testHomogeneousDirichlet( const std::string& meshFile, const uint_t& level, bool writeVTK = false )
{
   MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( mesh, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   P1DGEFunction< idx_t >         numerator( "numerator", storage, level, level );
   P1DGELaplaceOperator< real_t > L( storage, level, level );
   P1DGEMassOperator< real_t >    M( storage, level, level );

   numerator.enumerate( level );

   PETScSparseMatrix< P1DGELaplaceOperator< real_t > > Lpetsc;
   Lpetsc.createMatrixFromOperator( L, level, numerator, hyteg::All );

   PETScSparseMatrix< P1DGEMassOperator< real_t > > Mpetsc;
   Mpetsc.createMatrixFromOperator( M, level, numerator, hyteg::All );

   // solution as a lambda function
   std::function< real_t( const Point3D& p ) > u_x_expr = []( const Point3D& p ) -> real_t {
      const real_t x = p[0];
      const real_t y = p[1];
      return std::sin( M_PI * x ) * std::sin( M_PI * y ) * std::sin( M_PI * ( x + y ) );
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

   P1DGEFunction< real_t > u( "u", storage, level, level );
   P1DGEFunction< real_t > f( "f", storage, level, level );
   P1DGEFunction< real_t > rhs( "rhs", storage, level, level );
   P1DGEFunction< real_t > sol( "sol", storage, level, level );
   P1DGEFunction< real_t > err( "err", storage, level, level );
   P1DGEFunction< real_t > Merr( "Merr", storage, level, level );

   sol.interpolate( { u_x_expr, u_x_expr }, level, All );
   f.interpolate( { f_x_expr, f_x_expr }, level, All );

   // Why do I need this?
   f.interpolate( { u_x_expr, u_x_expr }, level, DirichletBoundary );

   // simulate a multiply operation:
   //    M.apply( f, rhs, level, All, Replace );
   auto f_vec = std::make_shared< PETScVector< real_t, P1DGEFunction > >();
   f_vec->createVectorFromFunction( f, numerator, level, All );
   auto rhs_vec = std::make_shared< PETScVector< real_t, P1DGEFunction > >();
   rhs_vec->createVectorFromFunction( f, numerator, level, All );
   Mpetsc.multiply( *f_vec, *rhs_vec );
   rhs_vec->createFunctionFromVector( rhs, numerator, level, All );

   rhs.interpolate(0, level, DirichletBoundary);
   u.getConformingPart()->interpolate( 0, level, DirichletBoundary );

   PETScCGSolver< P1DGELaplaceOperator< real_t > > solver( storage, level, numerator );
   solver.solve( L, u, rhs, level );

   err.assign( { 1.0, -1.0 }, { u, sol }, level );
   // calculate the error in the L2 norm
   auto err_vec = std::make_shared< PETScVector< real_t, P1DGEFunction > >();
   err_vec->createVectorFromFunction( err, numerator, level, All );
   auto Merr_vec = std::make_shared< PETScVector< real_t, P1DGEFunction > >();
   Merr_vec->createVectorFromFunction( err, numerator, level, All );
   Mpetsc.multiply( *err_vec, *Merr_vec );
   Merr_vec->createFunctionFromVector( Merr, numerator, level, All );

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

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();
   hyteg::PETScManager petscManager( &argc, &argv );

   for ( uint_t level = 2; level <= 8; level++ )
   {
      WALBERLA_LOG_INFO( hyteg::testHomogeneousDirichlet( "../../data/meshes/tri_1el.msh", level, true ) );
   }

   return 0;
}
