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
#include <cfenv>
#include <core/Environment.h>
#include <core/math/Constants.h>
#include <core/timing/Timer.h>

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

template < typename opType, template < class > class funcType >
void solveProblem( std::shared_ptr< hyteg::PrimitiveStorage >& storage, uint_t level, uint_t verbosity )
{
   FunctionTrait< funcType< real_t > > ft;
   std::string                         fileName;
   if ( ft.getTypeName().compare( "P2Function" ) == 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Running with P2 elements!" );
      fileName = "annulus_P2elements";
   }
   else if ( ft.getTypeName().compare( "P1Function / VertexDoFFunction" ) == 0 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Running with P1 elements!" );
      fileName = "annulus_P1elements";
   }
   else
   {
      WALBERLA_ABORT( "Houston we're in trouble here!" );
   }

   // Analytic solution in cartesian coordinates
   std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
      real_t m   = 5.0;
      real_t rho = std::sqrt( x[0] * x[0] + x[1] * x[1] );
      real_t phi = std::atan2( x[1], x[0] ) + pi;
      return std::pow( 2, m ) / ( std::pow( 2, 2 * m ) + 1 ) * ( std::pow( rho, m ) + std::pow( rho, -m ) ) * std::sin( m * phi );
   };
   funcType< real_t > u_exact( "u_analytic", storage, level, level );
   u_exact.interpolate( solFunc, level );

   // Create function for numeric solution and set Dirichlet boundary conditions
   funcType< real_t > u( "u_numeric", storage, level, level );
   u.interpolate( real_c( 0.0 ), level, Inner );
   u.interpolate( solFunc, level, DirichletBoundary );

   // Specify right-hand side of problem
   funcType< real_t > rhs( "rhs", storage, level, level );
   rhs.interpolate( real_c( 0.0 ), level, All );

   // Operator for weak-form
   opType lapOp( storage, level, level );

   // determine indices and dimensions
   funcType< idx_t > enumerator( "enumerator", storage, level, level );
   enumerator.enumerate( level );

   typedef typename FunctionTrait< funcType< idx_t > >::Tag    enumTag;
   uint_t                                                      globalDoFs = numberOfGlobalDoFs< enumTag >( *storage, level );
   uint_t                                                      localDoFs  = numberOfLocalDoFs< enumTag >( *storage, level );

   // assemble matrices
   PETScSparseMatrix< opType > lapMat( "Mat" );

   switch ( verbosity )
   {
   case 2: {
      lapMat.createMatrixFromOperator( lapOp, level, enumerator, All );
      lapMat.applyDirichletBC( enumerator, level );

      MatInfo info;
      MatGetInfo( lapMat.get(), MAT_GLOBAL_SUM, &info );
      WALBERLA_LOG_INFO_ON_ROOT( "Info on PETSc matrix:" );
      WALBERLA_LOG_INFO_ON_ROOT( "* block size ............................. " << real_c( info.block_size ) );
      WALBERLA_LOG_INFO_ON_ROOT( "* number of nonzeros ..................... " << info.nz_allocated );
      WALBERLA_LOG_INFO_ON_ROOT( "* memory allocated ....................... " << info.memory );
      WALBERLA_LOG_INFO_ON_ROOT( "* no. of matrix assemblies called ........ " << info.assemblies );
      WALBERLA_LOG_INFO_ON_ROOT( "* no. of mallocs during MatSetValues() ... " << info.mallocs << "\n" );

      // WALBERLA_LOG_INFO_ON_ROOT( "* no. of global DoFs (HyTeG) ............. " << globalDoFs );
      // WALBERLA_LOG_INFO_ON_ROOT( "* no. of local DoFs (HyTeG) .............. " << localDoFs );
      // uint_t globalDoFCount = numberOfGlobalInnerDoFs< enumTag >( *storage, level );
      // WALBERLA_LOG_INFO_ON_ROOT( "* no. of global inner DoFs (HyTeG) ....... "
      //                             << globalDoFCount );
      // break;
   }
      [[fallthrough]];
   case 1: {
      WALBERLA_LOG_INFO_ON_ROOT( "* no. of global DoFs (HyTeG) ............. " << globalDoFs );
      WALBERLA_LOG_INFO_ON_ROOT( "* no. of local DoFs (HyTeG) .............. " << localDoFs );
      uint_t globalDoFCount = numberOfGlobalInnerDoFs< enumTag >( *storage, level );
      WALBERLA_LOG_INFO_ON_ROOT( "* no. of global inner DoFs (HyTeG) ....... " << globalDoFCount );
      break;
   }
   default: {
   }
   }

   WALBERLA_LOG_INFO_ON_ROOT( " Entering PETSc solve " );

   // let PETSc solve the problem
   PETScLUSolver< opType > solver( storage, level );
   solver.solve( lapOp, u, rhs, level );

   // check size of error
   funcType< real_t > error( "error", storage, level, level );
   error.assign( { 1.0, -1.0 }, { u_exact, u }, level, All );

   VTKOutput vtkOutput( "../output", fileName, storage );
   // VTKOutput vtkOutput( "../output", "annulus", storage );
   vtkOutput.add( u_exact );
   vtkOutput.add( u );
   vtkOutput.add( error );
   vtkOutput.write( level );
}

int main( int argc, char* argv[] )
{
#ifndef __APPLE__
   #ifndef _MSC_VER
      feenableexcept( FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW );
   #endif
#endif
   // -------
   //  Setup
   // -------

   uint_t level     = 0;
   uint_t verbosity = 1;

   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // mesh, storage and geometry map for blending
   MeshInfo meshInfo = MeshInfo::meshAnnulus( 1.0, 2.0, 0.0, 2.0 * pi, MeshInfo::CROSS, 12, 2 );
   // MeshInfo meshInfo = MeshInfo::meshAnnulus( 1.0, 2.0, 0.0, 2.0 * pi, MeshInfo::CROSS, 6, 1 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   for ( auto it : setupStorage.getFaces() )
   {
      Face& face = *it.second;
      setupStorage.setGeometryMap( face.getID(), std::make_shared< AnnulusMap >( face ) );
   }
   for ( auto it : setupStorage.getEdges() )
   {
      Edge& edge = *it.second;
      setupStorage.setGeometryMap( edge.getID(), std::make_shared< AnnulusMap >( edge, setupStorage ) );
   }
   loadbalancing::roundRobin( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   if ( verbosity > 1 )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "" << setupStorage );
   }

   PETScManager petscManager( &argc, &argv );

   // --------------------------------
   //  Problem specification/solution
   // --------------------------------

   level = 5;
   solveProblem< P1ElementwiseBlendingLaplaceOperator, P1Function >( storage, level, verbosity );

   level = 4;
   solveProblem< P2ElementwiseBlendingLaplaceOperator, P2Function >( storage, level, verbosity );

   return 0;
}
