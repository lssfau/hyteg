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
real_t
    VectorLaplace( MeshInfo meshInfo, const uint_t& level, LambdaTuple sol_tuple, LambdaTuple rhs_tuple, bool writeVTK = false )
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
   }
   else
   {
      sol.interpolate( { u_x_expr, u_y_expr }, level, All );
      f.interpolate( { f_x_expr, f_y_expr }, level, All );
   }

   u.getConformingPart()->interpolate( 0, level, DirichletBoundary );
   M.apply( f, rhs, level, All, Replace );
   //rhs.interpolate( 0, level, DirichletBoundary );

   /*
   std::random_device                                         dev;
   std::mt19937                                               rng( dev() );
   std::uniform_int_distribution< std::mt19937::result_type > dist( 0, 1 );
   std::function< real_t( const Point3D& ) >                  randfunc = [&rng, &dist]( const Point3D& x ) { return dist( rng ); };
   u.interpolate( randfunc, level, Inner );
   */

   VTKOutput vtk( "../../output", "EGVectorLaplaceConvergenceTest", storage );
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
   if ( writeVTK )
   {
      vtk.write( level, 0 );
   }

   PETScCGSolver< EGLaplaceOperator > solver( storage, level, numerator );

   //CGSolver< EGLaplaceOperator > solver( storage, level, level );
   solver.solve( L, u, rhs, level );

   err.assign( { 1.0, -1.0 }, { u, sol }, level );

   // calculate the error in the L2 norm
   M.apply( err, Merr, level, All, Replace );
   auto discrL2 = sqrt( err.dotGlobal( Merr, level, Inner ) );
   //auto discrL2 = sqrt( err.dotGlobal( err, level, Inner ) / real_c( numberOfGlobalDoFs( u, level ) ) );

   if ( writeVTK )
   {
      vtk.write( level, 1 );
   }

   return discrL2;
}

void runTestcase( const std::string& name,
                  const uint_t&      minLevel,
                  const uint_t&      maxLevel,
                  MeshInfo           meshInfo,
                  LambdaTuple        sol_tuple,
                  LambdaTuple        rhs_tuple )
{
   auto l2ConvRate  = std::pow( 2, -( int( 1 ) + 1 ) );
   auto convRateEps = l2ConvRate * 0.1;
   auto err         = VectorLaplace( meshInfo, minLevel, sol_tuple, rhs_tuple );
   WALBERLA_LOG_INFO_ON_ROOT( "degree " << 1 << ", expected L2 rate: " << l2ConvRate
                                        << ", threshold: " << l2ConvRate + convRateEps );
   WALBERLA_LOG_INFO_ON_ROOT( "error level " << minLevel << ": " << err );
   for ( uint_t l = minLevel + 1; l <= maxLevel; l++ )
   {
      auto errFiner     = VectorLaplace( meshInfo, l, sol_tuple, rhs_tuple, true );
      auto computedRate = errFiner / err;

      WALBERLA_LOG_INFO_ON_ROOT( "error level " << l << ": " << errFiner );
      WALBERLA_LOG_INFO_ON_ROOT( "computed rate level " << l << " / " << l - 1 << ": " << computedRate );

      /* WALBERLA_CHECK_LESS_EQUAL( computedRate,
                                 l2ConvRate + convRateEps,
                                 "Convergence L2 rate level " << l << " vs level " << l - 1
                                                              << " not sufficiently small (computed: " << computedRate
                                                              << ", estimated + eps: " << l2ConvRate + convRateEps << ")" );*/
      WALBERLA_LOG_INFO_ON_ROOT( "Convergence L2 rate level " << l << " vs level " << l - 1
                                                              << " not sufficiently small (computed: " << computedRate
                                                              << ", estimated + eps: " << l2ConvRate + convRateEps << ")" );

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

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single face, hom. BC, rhs != 0 ###" );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
         return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * ( x[0] + x[1] - 1 ) );
      };

      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
         return 4 * pi * pi * ( -2 * sin( 4 * pi * ( x[0] + x[1] ) ) + sin( 4 * pi * x[0] ) + sin( 4 * pi * x[1] ) );
      };
      hyteg::runTestcase( "1face_Dirichlet0",
                          3,
                          6,
                          MeshInfo::meshFaceChain( 1 ),
                          std::make_tuple( solFunc, solFunc, solFunc ),
                          std::make_tuple( rhsFunc, rhsFunc, rhsFunc ) );
   }

  

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single tetrahedron, hom. BC, rhs != 0 ###" );

      std::function< real_t( const Point3D& ) > one_tet_sol = []( const Point3D& x ) {
         return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] ) * sin( 2 * pi * ( x[0] + x[1] + x[2] - 1 ) );
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
                          3,
                          6,
                          MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" ),
                          std::make_tuple( one_tet_sol, one_tet_sol, one_tet_sol ),
                          std::make_tuple( one_tet_rhs, one_tet_rhs, one_tet_rhs ) );
   }
   /*
 {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single triangle, inhom. BC, rhs = 0 ###" );

      MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& ) { return 0; };

      hyteg::runTestcase( "1tri",
                          3,
                          6,
                          MeshInfo::fromGmshFile( "../../data/meshes/tri_1el.msh" ),
                          std::make_tuple( solFunc, solFunc, solFunc ),
                          std::make_tuple( rhsFunc, rhsFunc, rhsFunc ) );
   }
   */
   /*
      hyteg::run( "Cube_Dirichlet0",
      3,7,
                  MeshInfo::meshCuboid( Point3D( { -1, -1, -1 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 ),
                  std::make_tuple(
                      []( const Point3D& p ) -> real_t {
                         const real_t x = p[0];
                         const real_t y = p[1];
                         const real_t z = p[2];
                         return std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) *
                                std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) *
                                std::sin( M_PI * ( ( 1.0 / 2.0 ) * z + 1.0 / 2.0 ) );
                      },
                      []( const Point3D& p ) -> real_t {
                         const real_t x = p[0];
                         const real_t y = p[1];
                         const real_t z = p[2];
                         return std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) *
                                std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) *
                                std::sin( M_PI * ( ( 1.0 / 2.0 ) * z + 1.0 / 2.0 ) );
                      },
                      []( const Point3D& p ) -> real_t {
                         const real_t x = p[0];
                         const real_t y = p[1];
                         const real_t z = p[2];
                         return std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) *
                                std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) *
                                std::sin( M_PI * ( ( 1.0 / 2.0 ) * z + 1.0 / 2.0 ) );
                      } ),
                  std::make_tuple(
                      []( const Point3D& p ) -> real_t {
                         const real_t x = p[0];
                         const real_t y = p[1];
                         const real_t z = p[2];
                         return -3.0 / 4.0 * std::pow( M_PI, 2 ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) *
                                std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) *
                                std::sin( M_PI * ( ( 1.0 / 2.0 ) * z + 1.0 / 2.0 ) );
                      },
                      []( const Point3D& p ) -> real_t {
                         const real_t x = p[0];
                         const real_t y = p[1];
                         const real_t z = p[2];
                         return -3.0 / 4.0 * std::pow( M_PI, 2 ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) *
                                std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) *
                                std::sin( M_PI * ( ( 1.0 / 2.0 ) * z + 1.0 / 2.0 ) );
                      },
                      []( const Point3D& p ) -> real_t {
                         const real_t x = p[0];
                         const real_t y = p[1];
                         const real_t z = p[2];
                         return -3.0 / 4.0 * std::pow( M_PI, 2 ) * std::sin( M_PI * ( ( 1.0 / 2.0 ) * x + 1.0 / 2.0 ) ) *
                                std::sin( M_PI * ( ( 1.0 / 2.0 ) * y + 1.0 / 2.0 ) ) *
                                std::sin( M_PI * ( ( 1.0 / 2.0 ) * z + 1.0 / 2.0 ) );
                      } ) );
   }
   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on single macro, hom. BC, rhs != 0 ###" );

      MeshInfo meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/3D/tet_1el.msh" );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
         return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] ) * sin( 2 * pi * ( x[0] + x[1] + x[2] - 1 ) );
      };

      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
         return -( -24 * pi * pi * sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] ) *
                       sin( 2 * pi * ( x[0] + x[1] + x[2] - 1 ) ) +
                   8 * pi * pi * sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * cos( 2 * pi * x[2] ) *
                       cos( 2 * pi * ( x[0] + x[1] + x[2] - 1 ) ) +
                   8 * pi * pi * cos( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] ) *
                       cos( 2 * pi * ( x[0] + x[1] + x[2] - 1 ) ) +
                   8 * pi * pi * sin( 2 * pi * x[0] ) * cos( 2 * pi * x[1] ) * sin( 2 * pi * x[2] ) *
                       cos( 2 * pi * ( x[0] + x[1] + x[2] - 1 ) ) );
      };

      hyteg::runTest( 3, 5, 1, 1, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, hom. BC, rhs != 0 ###" );

      MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
         return sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] );
      };

      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
         return 12 * pi * pi * ( sin( 2 * pi * x[0] ) * sin( 2 * pi * x[1] ) * sin( 2 * pi * x[2] ) );
      };

      hyteg::runTest( 2, 4, 1, 1, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, inhom. BC, rhs = 0 ###" );

      MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) { return sin( x[0] ) * sinh( x[1] ) * x[2]; };
      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& ) { return 0; };

      hyteg::runTest( 2, 4, 1, 1, meshInfo, solFunc, rhsFunc );
   }

   WALBERLA_LOG_INFO_ON_ROOT( "" );

   {
      WALBERLA_LOG_INFO_ON_ROOT( "### Test on multiple macros, inhom. BC, rhs != 0 ###" );

      MeshInfo meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 1, 1, 1 );

      std::function< real_t( const Point3D& ) > solFunc = []( const Point3D& x ) {
         return sin( x[0] ) * sin( x[1] ) * sin( x[2] );
      };

      std::function< real_t( const Point3D& ) > rhsFunc = []( const Point3D& x ) {
         return 3 * sin( x[0] ) * sin( x[1] ) * sin( x[2] );
      };

      hyteg::runTest( 2, 4, 1, 1, meshInfo, solFunc, rhsFunc );
   }*/

   return EXIT_SUCCESS;
}
