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

#include <random>

#include "core/DataTypes.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P1DGEP0StokesFunction.hpp"
#include "hyteg/composites/P1DGEP0StokesOperator.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodStokesOperator.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/egfunctionspace/EGOperators.hpp"
#include "hyteg/functions/FunctionTraits.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.cpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/petsc/PETScMinResSolver.hpp"
#include "hyteg/petsc/PETScSparseMatrix.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using walberla::real_t;
using walberla::uint_t;

using hyteg::dg::eg::EGLaplaceOperator;
using hyteg::dg::eg::EGMassOperator;
using hyteg::dg::eg::EGP0ConstEpsilonStokesOperator;
using hyteg::dg::eg::EGP0StokesOperator;

// scalar lambda for one component of analytical solution and rhs
typedef std::function< real_t( const hyteg::PointND< real_t, 3 >& p ) > ScalarLambda;

// tuple of function for solution (u,p) and rhs of vector values stokes equation
typedef std::tuple< ScalarLambda, ScalarLambda, ScalarLambda > LambdaTuple;

namespace hyteg {
real_t VectorLaplace( const std::string& name,
                      MeshInfo           meshInfo,
                      const uint_t&      level,
                      LambdaTuple        sol_tuple,
                      LambdaTuple        rhs_tuple,
                      int                solverType,
                      bool               writeVTK = false )
{
   // MeshInfo              mesh = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   EGFunction< idx_t > numerator( "numerator", storage, level, level );
   numerator.enumerate( level );
   EGLaplaceOperator L( storage, level, level );
   EGMassOperator    M( storage, level, level );

   EGFunction< real_t > u( "u", storage, level, level );
   EGFunction< real_t > f( "f", storage, level, level );
   EGFunction< real_t > rhs( "rhs", storage, level, level );
   EGFunction< real_t > sol( "sol", storage, level, level );
   EGFunction< real_t > err( "err", storage, level, level );
   EGFunction< real_t > Merr( "Merr", storage, level, level );

   WALBERLA_LOG_INFO_ON_ROOT( "dofs: " << u.getNumberOfGlobalDoFs( level ) );

   auto [u_x_expr, u_y_expr, u_z_expr] = sol_tuple;
   auto [f_x_expr, f_y_expr, f_z_expr] = rhs_tuple;
   if ( storage->hasGlobalCells() )
   {
      sol.interpolate( { u_x_expr, u_y_expr, u_z_expr }, level, All );
      f.interpolate( { f_x_expr, f_y_expr, f_z_expr }, level, All );
      u.getConformingPart()->interpolate( { u_x_expr, u_y_expr, u_z_expr }, level, DirichletBoundary );
   }
   else
   {
      sol.interpolate( { u_x_expr, u_y_expr }, level, All );
      f.interpolate( { f_x_expr, f_y_expr }, level, All );
      u.getConformingPart()->interpolate( { u_x_expr, u_y_expr }, level, DirichletBoundary );
   }

   M.apply( f, rhs, level, Inner, Replace );
   rhs.interpolate( 0, level, DirichletBoundary );

   /*
   std::random_device                                         dev;
   std::mt19937                                               rng( dev() );
   std::uniform_int_distribution< std::mt19937::result_type > dist( 0, 1 );
   std::function< real_t( const Point3D& ) >                  randfunc = [&rng, &dist]( const Point3D& x ) { return dist( rng ); };
   u.interpolate( randfunc, level, Inner );
   */

   switch ( solverType )
   {
   case 0: {
      PETScCGSolver< EGLaplaceOperator > solver( storage, level, numerator );
      solver.solve( L, u, rhs, level );
      break;
   }
   case 1:{
      CGSolver< EGLaplaceOperator > solver( storage, level, level );
      solver.solve( L, u, rhs, level );
      break;}
   default:
      WALBERLA_ABORT("No solver chosen.");
   }

   // u.getDiscontinuousPart()->interpolate( 0, level, DirichletBoundary );
   err.assign( { 1.0, -1.0 }, { u, sol }, level );

   // calculate the error in the L2 norm
   M.apply( err, Merr, level, All, Replace );
   auto discrL2 = sqrt( err.dotGlobal( Merr, level, Inner ) );
   //auto discrL2 = sqrt( err.dotGlobal( err, level, Inner ) / real_c( numberOfGlobalDoFs( u, level ) ) );

   if ( writeVTK )
   {
      VTKOutput vtk( "../../output", name, storage );
      vtk.add( u );
      vtk.add( *u.getConformingPart() );
      vtk.add( *u.getDiscontinuousPart() );
      vtk.add( sol );
      vtk.add( *sol.getConformingPart() );
      vtk.add( *sol.getDiscontinuousPart() );
      vtk.add( f );
      vtk.add( *f.getConformingPart() );
      vtk.add( *f.getDiscontinuousPart() );
      vtk.add( rhs );
      vtk.add( *rhs.getConformingPart() );
      vtk.add( *rhs.getDiscontinuousPart() );

      vtk.add( err );
      vtk.add( *err.getConformingPart() );
      vtk.add( *err.getDiscontinuousPart() );
      vtk.write( level );
   }

   return discrL2;
}

void runTestcase( const std::string& name,
                  const uint_t&      minLevel,
                  const uint_t&      maxLevel,
                  MeshInfo           meshInfo,
                  LambdaTuple        sol_tuple,
                  LambdaTuple        rhs_tuple,
                  int                solverType )
{
   auto l2ConvRate  = std::pow( 2, -( int( 1 ) + 1 ) );
   auto convRateEps = l2ConvRate * 0.1;
   auto err         = VectorLaplace( name, meshInfo, minLevel, sol_tuple, rhs_tuple, solverType );
   WALBERLA_LOG_INFO_ON_ROOT( "degree " << 1 << ", expected L2 rate: " << l2ConvRate
                                        << ", threshold: " << l2ConvRate + convRateEps );
   WALBERLA_LOG_INFO_ON_ROOT( "error level " << minLevel << ": " << err );
   for ( uint_t l = minLevel + 1; l <= maxLevel; l++ )
   {
      auto errFiner     = VectorLaplace( name, meshInfo, l, sol_tuple, rhs_tuple, solverType, true );
      auto computedRate = errFiner / err;

      WALBERLA_LOG_INFO_ON_ROOT( "error level " << l << ": " << errFiner );
      WALBERLA_LOG_INFO_ON_ROOT( "computed rate level " << l << " / " << l - 1 << ": " << computedRate );

      /* WALBERLA_CHECK_LESS_EQUAL( computedRate,
                                 l2ConvRate + convRateEps,
                                 "Convergence L2 rate level " << l << " vs level " << l - 1
                                                              << " not sufficiently small (computed: " << computedRate
                                                              << ", estimated + eps: " << l2ConvRate + convRateEps << ")" );*/
      WALBERLA_LOG_INFO_ON_ROOT( "Convergence L2 rate level " << l << " vs level " << l - 1 << "; computed: " << computedRate
                                                              << ", estimated + eps: " << l2ConvRate + convRateEps );

      err = errFiner;
   }
}
} // namespace hyteg

int main( int argc, char** argv )
{
   using hyteg::MeshInfo;
   using hyteg::Point3D;
   using walberla::real_t;

   using walberla::math::pi;
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PETScManager petscManager( &argc, &argv );
   int                 minLevel = atoi( argv[3] );
   int                 maxLevel = atoi( argv[4] );
   int                 solverType   = atoi( argv[5] );
   switch ( atoi( argv[2] ) )
   {
   case 0:
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single face, hom. BC, rhs != 0 ###" );
      {
         std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
            return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * ( x[0] + x[1] - 1 ) );
         };

         std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
            return 4 * pi * pi * ( -2 * sin( 4 * pi * ( x[0] + x[1] ) ) + sin( 4 * pi * x[0] ) + sin( 4 * pi * x[1] ) );
         };
         hyteg::runTestcase( "1tri_Dirichlet0",
                             minLevel,
                             maxLevel,
                             MeshInfo::meshFaceChain( 1 ),
                             std::make_tuple( solFunc, solFunc, solFunc ),
                             std::make_tuple( rhsFunc, rhsFunc, rhsFunc ),
                             solverType );
      }
      break;
   case 1:
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single triangle, inhom. BC, rhs = 0 ###" );
      {
         std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
         std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& ) { return 0; };

         hyteg::runTestcase( "1tri",
                             minLevel,
                             maxLevel,
                             MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ),
                             std::make_tuple( solFunc, solFunc, solFunc ),
                             std::make_tuple( rhsFunc, rhsFunc, rhsFunc ),
                             solverType );
      }
      break;
   case 2:
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple triangles, inhom. BC, rhs = 0 ###" );
      {
         MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/quad_4el.msh" );

         std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
         std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& ) { return 0; };

         hyteg::runTestcase( "quad_4el",
                             minLevel,
                             maxLevel,
                             meshInfo,
                             std::make_tuple( solFunc, solFunc, solFunc ),
                             std::make_tuple( rhsFunc, rhsFunc, rhsFunc ),
                             solverType );
      }
      break;
   case 3:
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on one tet, hom. BC, rhs != 0 ###" );
      {
         std::function< real_t( const Point3D& ) > one_tet_sol = []( const Point3D& x ) {
            return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] ) *
                   sin( 2 * pi * ( x[0] + x[1] + x[2] - 1 ) );
         };

         std::function< real_t( const Point3D& ) > one_tet_rhs = []( const Point3D& x ) {
            return -( -24 * pi * pi * sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] ) *
                          sin( 2 * pi * ( x[0] + x[1] + x[2] - 1 ) ) +
                      8 * pi * pi * sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * cos( 2 * pi * x[2] ) *
                          cos( 2 * pi * ( x[0] + x[1] + x[2] - 1 ) ) +
                      8 * pi * pi * cos( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] ) *
                          cos( 2 * pi * ( x[0] + x[1] + x[2] - 1 ) ) +
                      8 * pi * pi * sin( 2 * pi * x[0] ) * cos( 2 * pi * x[1] ) * sin( 2 * pi * x[2] ) *
                          cos( 2 * pi * ( x[0] + x[1] + x[2] - 1 ) ) );
         };

         hyteg::runTestcase( "1tet_Dirichlet0",
                             minLevel,
                             maxLevel,
                             MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ),
                             std::make_tuple( one_tet_sol, one_tet_sol, one_tet_sol ),
                             std::make_tuple( one_tet_rhs, one_tet_rhs, one_tet_rhs ),
                             solverType );
      }
      break;
   case 4:
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on one tet, inhom. BC, rhs != 0 ###" );
      {
         std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
            return sin( x[0] ) * sin( x[1] ) * sin( x[2] );
         };

         std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
            return 3 * sin( x[0] ) * sin( x[1] ) * sin( x[2] );
         };

         hyteg::runTestcase( "1tet",
                             minLevel,
                             maxLevel,
                             MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ),
                             std::make_tuple( solFunc, solFunc, solFunc ),
                             std::make_tuple( rhsFunc, rhsFunc, rhsFunc ),
                             solverType );
      }
      break;
   case 5:
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on one tet, inhom. BC, rhs != 0 (second) ###" );
      {
         std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
            return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] );
         };

         std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
            return 12 * pi * pi * ( sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] ) );
         };

         hyteg::runTestcase( "1tet_second",
                             minLevel,
                             maxLevel,
                             MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ),
                             std::make_tuple( solFunc, solFunc, solFunc ),
                             std::make_tuple( rhsFunc, rhsFunc, rhsFunc ),
                             solverType );
      }
      break;
   case 6:
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on cube, hom. BC, rhs != 0 ###" );
      {
         MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );

         std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
            return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] );
         };

         std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
            return 12 * pi * pi * ( sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] ) );
         };

         hyteg::runTestcase( "cube_Dirichlet0",
                             minLevel,
                             maxLevel,
                             meshInfo,
                             std::make_tuple( solFunc, solFunc, solFunc ),
                             std::make_tuple( rhsFunc, rhsFunc, rhsFunc ),
                             solverType );
      }
      break;
   case 7:
      WALBERLA_LOG_INFO_ON_ROOT( "### Test nonzero on interfaces, on cube, hom. BC, rhs != 0 ###" );
      {
         MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );

         std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
            return sin( pi * x[0] ) * sin( pi * x[1] ) * sin( pi * x[2] );
         };

         std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
            return 3 * pi * pi * ( sin( pi * x[0] ) * sin( pi * x[1] ) * sin( pi * x[2] ) );
         };

         hyteg::runTestcase( "cube_Dirichlet0_nonzero_on_interfaces",
                             minLevel,
                             maxLevel,
                             meshInfo,
                             std::make_tuple( solFunc, solFunc, solFunc ),
                             std::make_tuple( rhsFunc, rhsFunc, rhsFunc ),
                             solverType );
      }
      break;
   case 8:
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on cube, inhom. BC, rhs = 0 ###" );
      {
         MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );

         std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) { return sin( x[0] ) * sinh( x[1] ) * x[2]; };
         std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& ) { return 0; };
         hyteg::runTestcase( "cube_rhs0",
                             minLevel,
                             maxLevel,
                             meshInfo,
                             std::make_tuple( solFunc, solFunc, solFunc ),
                             std::make_tuple( rhsFunc, rhsFunc, rhsFunc ),
                             solverType );
      }
      break;
   case 9:
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on cube, inhom. BC, rhs != 0 second ###" );
      {
         MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );

         std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
            return sin( x[0] ) * sin( x[1] ) * sin( x[2] );
         };

         std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
            return 3 * sin( x[0] ) * sin( x[1] ) * sin( x[2] );
         };

         hyteg::runTestcase( "cube",
                             minLevel,
                             maxLevel,
                             meshInfo,
                             std::make_tuple( solFunc, solFunc, solFunc ),
                             std::make_tuple( rhsFunc, rhsFunc, rhsFunc ),
                             solverType );
      }
      break;
   case 10:
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on pyramid of 2 elements, inhom. BC, rhs != 0 ###" );
      {
         MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/pyramid_2el.msh" );

         std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
            return sin( x[0] ) * sin( x[1] ) * sin( x[2] );
         };

         std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
            return 3 * sin( x[0] ) * sin( x[1] ) * sin( x[2] );
         };

         hyteg::runTestcase( "pyramid_2el",
                             minLevel,
                             maxLevel,
                             meshInfo,
                             std::make_tuple( solFunc, solFunc, solFunc ),
                             std::make_tuple( rhsFunc, rhsFunc, rhsFunc ),
                             solverType );
      }
      break;
   default:
      WALBERLA_ABORT( "No testcase chosen" );
   }

   return EXIT_SUCCESS;
}
