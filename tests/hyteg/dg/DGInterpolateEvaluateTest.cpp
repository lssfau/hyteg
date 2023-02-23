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

#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/dgfunctionspace/DGMassForm_Example.hpp"
#include "hyteg/dgfunctionspace/DGOperator.hpp"
#include "hyteg/petsc/PETScCGSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::math::pi;

void test( int dim, uint_t level, uint_t degree, std::function< real_t( const Point3D& ) > f, real_t maxPointwiseError )
{
   using namespace dg;

   walberla::math::seedRandomGenerator( 12345678 );
   const uint_t numRandomEvaluations = 1000;

   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( { 0, 0 } ), Point2D( { 1, 1 } ), MeshInfo::CRISS, 2, 2 );
   if ( dim == 3 )
   {
      meshInfo = MeshInfo::meshCuboid( Point3D( { 0, 0, 0 } ), Point3D( { 1, 1, 1 } ), 2, 2, 2 );
   }

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, 1 );

   auto basis    = std::make_shared< DGBasisLinearLagrange_Example >();
   auto massForm = std::make_shared< DGMassForm_Example >();

   P1Function< real_t > interpolatedP1( "u_P1", storage, level, level );
   DGFunction< real_t > u( "u", storage, level, level, basis, degree );
   DGFunction< real_t > tmp( "tmp", storage, level, level, basis, degree );

   DGFunction< idx_t > numerator( "numerator", storage, level, level, basis, degree );
   numerator.enumerate( level );

   DGOperator M( storage, level, level, massForm );

   // Interpolate solution into u
   tmp.evaluateLinearFunctional( f, level );
   PETScCGSolver< DGOperator > solverM( storage, level, numerator );
   solverM.solve( M, u, tmp, level );

   // Interpolate solution for P1
   interpolatedP1.interpolate( f, level, All );

   VTKOutput vtk( "../../output/", "DGInterpolateEvaluateTest", storage );
   vtk.add( u );
   vtk.add( tmp );
   vtk.add( interpolatedP1 );
   vtk.write( level );

   for ( uint_t i = 0; i < numRandomEvaluations; ++i )
   {
      Point3D coordinates(Point3D::Zero());
      coordinates[0] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );
      coordinates[1] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );
      coordinates[2] = real_c( walberla::math::realRandom( 0.0, 1.0 ) );

      real_t     value;
      const bool success = u.evaluate( coordinates, level, value );
      WALBERLA_CHECK( success, "Could not evaluate successfully." );

      const real_t err = std::abs( value - f( coordinates ) );

      WALBERLA_CHECK_LESS(
          err, maxPointwiseError, "Failed at " << coordinates << ", evaluated value: " << value << ", error: " << err );
   }
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_CHECK_EQUAL(
       walberla::mpi::MPIManager::instance()->numProcesses(), 1, "Evaluate test only implemented for a single process." );

   hyteg::PETScManager petscManager( &argc, &argv );

   hyteg::test(
       2, 4, 1, []( const hyteg::Point3D& ) { return 1; }, 1e-11 );
   hyteg::test(
       2, 4, 1, []( const hyteg::Point3D& x ) { return x[0] - 2 * x[1]; }, 1e-11 );

   hyteg::test(
       3, 3, 1, []( const hyteg::Point3D& ) { return 1; }, 1e-11 );
   hyteg::test(
       3, 3, 1, []( const hyteg::Point3D& x ) { return x[0] - 2 * x[1] + 3 * x[2]; }, 1e-11 );

   return EXIT_SUCCESS;
}
