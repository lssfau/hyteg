/*
 * Copyright (c) 2023 Marcus Mohr.
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

// Check that application of Dirichlet Boundary conditions with the
// PETScLUSolver class works when using user-defined mesh boundary
// flags and BoundaryUID objects.

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

#ifndef HYTEG_BUILD_WITH_PETSC
#error "This test only works with PETSc enabled. Please enable it via -DHYTEG_BUILD_WITH_PETSC=ON"
#endif

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   PETScManager petscManager( &argc, &argv );

   const uint_t nx       = 1;
   const uint_t ny       = 1;
   const uint_t minLevel = 2;
   const uint_t maxLevel = 2;

   // generate mesh and setup storage
   auto meshInfo     = MeshInfo::meshRectangle( Point2D( 0, 0 ), Point2D( 1, 1 ), MeshInfo::CROSS, nx, ny );
   auto setupStorage = std::make_shared< SetupPrimitiveStorage >(
       meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   // callbacks for marking parts of the boundary
   real_t eps         = real_c( 1e-3 );
   auto   topBoundary = [eps]( const Point3D& x ) { return std::abs( x[1] - real_c( 1 ) ) < eps; };

   auto bottomBoundary = [eps]( const Point3D& x ) { return x[1] < eps; };

   auto leftBoundary = [eps]( const Point3D& x ) {
      return x[0] < eps && !( std::abs( x[1] - real_c( 1 ) ) < eps ) && !( x[1] < eps );
   };

   auto rightBoundary = [eps]( const Point3D& x ) {
      return std::abs( x[0] - real_c( 1 ) ) < eps && !( std::abs( x[1] - real_c( 1 ) ) < eps ) && !( x[1] < eps );
   };

   // assigning mesh boundary flags to different parts of the boundary
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( 1, topBoundary );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( 2, bottomBoundary );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( 3, leftBoundary );
   setupStorage->setMeshBoundaryFlagsByCentroidLocation( 4, rightBoundary );

   // now create the actual storage
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   // setup bounday condition information
   BoundaryCondition bc;

   BoundaryUID idTopWall    = bc.createDirichletBC( "topWall", { 1 } );
   BoundaryUID idBottomWall = bc.createDirichletBC( "botWall", { 2 } );
   BoundaryUID idLeftWall   = bc.createDirichletBC( "leftWall", { 3 } );
   BoundaryUID idRightWall  = bc.createDirichletBC( "rightWall", { 4 } );

   // using a quadratic polynomial
   std::function< real_t( const Point3D& ) > poly2 = []( const Point3D& x ) {
      return real_c( 2 ) * x[0] * x[0] + x[0] * x[1] + real_c( 0.5 ) * x[1] * x[1];
   };
   std::function< real_t( const Point3D& ) > poly2laplacian = []( const Point3D& ) { return real_c( -5 ); };

   // setup functions
   P1Function< real_t > uAnalytic( "analytical solution", storage, minLevel, maxLevel, bc );
   P1Function< real_t > uDiscrete( "discrete solution", storage, minLevel, maxLevel, bc );
   P1Function< real_t > difference( "difference", storage, minLevel, maxLevel, bc );
   P1Function< real_t > rhs( "rhs", storage, minLevel, maxLevel, bc );
   P1Function< real_t > weakRhs( "weak rhs", storage, minLevel, maxLevel, bc );

   uAnalytic.interpolate( poly2, maxLevel, All );
   uDiscrete.interpolate( poly2, maxLevel, idTopWall );
   uDiscrete.interpolate( poly2, maxLevel, idBottomWall );
   uDiscrete.interpolate( poly2, maxLevel, idLeftWall );
   uDiscrete.interpolate( poly2, maxLevel, idRightWall );

   // prepare right-hand side of weak formulation
   P1ConstantMassOperator massOperator( storage, minLevel, maxLevel );
   rhs.interpolate( poly2laplacian, maxLevel, All );
   massOperator.apply( rhs, weakRhs, maxLevel, All );

   // solve the linear system for the negative Laplacian
   P1ConstantLaplaceOperator                  diffusionOperator( storage, minLevel, maxLevel );
   PETScLUSolver< P1ConstantLaplaceOperator > solver( storage, maxLevel );
   solver.solve( diffusionOperator, uDiscrete, weakRhs, maxLevel );

   // compare results
   difference.assign( { real_c( 1 ), real_c( -1 ) }, { uAnalytic, uDiscrete }, maxLevel, All );
   real_t diffMax = difference.getMaxDoFMagnitude( maxLevel, All );

   // export functions to visually check results
   bool beVerbose = false;
   if ( beVerbose )
   {
      std::string fPath = "../../output";
      std::string fName = "PetscApplyDirichletBoundaryUIDTest";
      WALBERLA_LOG_INFO_ON_ROOT( "Exporting to '" << fPath << "/" << fName << "'" );
      VTKOutput vtkOutput( fPath, fName, storage );
      vtkOutput.add( uAnalytic );
      vtkOutput.add( uDiscrete );
      vtkOutput.add( difference );
      vtkOutput.write( maxLevel );
   }

   // perform the test
   real_t threshold = std::is_same< real_t, double >::value ? 1e-15 : 1e-7;
   WALBERLA_CHECK_LESS( diffMax, threshold );
}
