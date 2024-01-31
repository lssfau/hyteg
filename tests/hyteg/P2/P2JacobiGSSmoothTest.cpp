/*
 * Copyright (c) 2017-2019 Marcus Mohr.
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
#include <hyteg/communication/Syncing.hpp>
#include <hyteg/p1functionspace/VertexDoFMacroEdge.hpp>
#include <hyteg/p1functionspace/VertexDoFMacroFace.hpp>

#include "core/Environment.h"
#include "core/math/Random.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "constantStencilOperator/P2ConstantOperator.hpp"

using walberla::real_t;
using namespace hyteg;

// ----------------------------------------------------------------
//  Perform a single undamped Jacobi smoothing step on a quadratic
//  polynomial and check that this does not alter the polynomial;
//  do the same for Gauss-Seidel
// ----------------------------------------------------------------

// Coefficients of the test polynomial
#define COEFF_A00 real_c(  2.0 )
#define COEFF_A10 real_c(  0.0 )
#define COEFF_A01 real_c(  0.0 )
#define COEFF_A11 real_c(  1.0 )
#define COEFF_A20 real_c(  0.5 )
#define COEFF_A02 real_c( -1.0 )

// Do we want output for visualisation?
bool outputVTK = true;

// Auxilliary type for code de-duplication
typedef enum { JACOBI, GAUSS_SEIDEL } SmootherType;

// Jacobi smoothing
template<typename OperatorType>
void smooth( OperatorType& oper,
             P2Function< real_t >& uNew, 
             P2Function< real_t >& rhs, 
             P2Function< real_t >& uOld,
             uint_t level ) {
  oper.smooth_jac( uNew, rhs, uOld, 1.0, level, Inner );
}

// Gauss-Seidel smoothing
template<typename OperatorType>
void smooth( OperatorType& oper,
             P2Function< real_t >& u, 
             P2Function< real_t >& rhs, 
             uint_t level ) {
  oper.smooth_gs( u, rhs, level, Inner );
}

void smooth( P2ElementwiseLaplaceOperator& oper,
             P2Function< real_t >& u, 
             P2Function< real_t >& rhs, 
             uint_t level ) {
  WALBERLA_UNUSED( oper  );
  WALBERLA_UNUSED( u     );
  WALBERLA_UNUSED( rhs   );
  WALBERLA_UNUSED( level );
  WALBERLA_ABORT( "Elementwise operator does not offer Gauss-Seidel smoothing!" );
}


template<typename OperatorTypeMass, typename OperatorTypeLaplace >
void testSmoother( SmootherType type,
                   P2Function< real_t >& rhsStrong,
                   P2Function< real_t >& polynomial,
                   std::shared_ptr<PrimitiveStorage> storage,
                   uint_t level ) {

  // Setup operators
  OperatorTypeMass massOp( storage, level, level );
  OperatorTypeLaplace Lap( storage, level, level );

  // Determine right-hand side
  P2Function< real_t > rhs( "right-hand side", storage, level, level );
  massOp.apply( rhsStrong, rhs, level, Inner );
  rhs.interpolate( 0.0, level, DirichletBoundary );

  // Execute a single smoothing step
  P2Function< real_t > smoothed( "after one smoothing step", storage, level, level );
  switch( type ) {
  case JACOBI:
    smooth( Lap, smoothed, rhs, polynomial, level );
    break;
  case GAUSS_SEIDEL:
    smoothed.assign( {1.0}, {polynomial}, level, All );
    smooth( Lap, smoothed, rhs, level );
    break;
  }

  // Compute difference to before
  P2Function< real_t > difference( "difference", storage, level, level );
  difference.assign( {1.0, -1.0}, {polynomial, smoothed}, level, Inner );

  // output data for visualisation
  if ( outputVTK )
    {
      std::string filename;
      switch( type ) {
      case JACOBI:
        filename.assign( "P2_Jacobi_Smooth_Test" );
        break;
      case GAUSS_SEIDEL:
        filename.assign( "P2_GS_Smooth_Test" );
        break;
      }
      VTKOutput vtkOutput( "../../output", filename, storage );
      vtkOutput.add( difference );
      vtkOutput.add( smoothed   );
      vtkOutput.add( polynomial );
      vtkOutput.add( rhs        );
      vtkOutput.write( level );
   }

  // Check size of difference for Jacobi smoother
  uint_t nDoFs = numberOfGlobalDoFs<P2FunctionTag>( (*storage), level );
  real_t errorNorm = std::sqrt( difference.dotGlobal( difference, level ) / real_c(nDoFs) );
  auto dp = std::is_same< real_t, double >();
  WALBERLA_CHECK_LESS( errorNorm, dp ? 2e-13 : 1.5e-7 );

  errorNorm = difference.getMaxMagnitude( level );
  WALBERLA_CHECK_LESS( errorNorm, dp ? 3e-13 : 8e-7 );

  switch( type ) {
  case JACOBI:
    WALBERLA_LOG_INFO_ON_ROOT( "  Check passed for Jacobi" );
    break;
  case GAUSS_SEIDEL:
    WALBERLA_LOG_INFO_ON_ROOT( "  Check passed for Gauss-Seidel" );
    break;
  }
}


int main( int argc, char** argv )
{
   uint_t level     = 2;

   // Setup enviroment
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   // Generate a regular mesh of a square centered around the origin with a criss pattern and 18 triangles
   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( -1.0, -1.0 ), Point2D( 1.0, 1.0 ), MeshInfo::CRISS, 3, 3 );

   // Prepare primitive storage
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   loadbalancing::roundRobin( setupStorage );
   WALBERLA_LOG_INFO_ON_ROOT( setupStorage );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   // Setup polynomial
   std::function< real_t( const Point3D& ) > myPolyFunc = []( const Point3D& x ) {
      real_t val = COEFF_A00 + COEFF_A10 * x[0] + COEFF_A01 * x[1] + COEFF_A11 * x[0] * x[1];
      val += COEFF_A20 * x[0] * x[0] + COEFF_A02 * x[1] * x[1];
      return val;
   };
   P2Function< real_t > poly( "quadratic polynomial", storage, level, level );
   poly.interpolate( myPolyFunc, level );

   // Setup strong right-hand side
   P2ConstantMassOperator massOp( storage, level, level );
   P2Function< real_t >   rhsStrong( "right-hand side", storage, level, level );
   rhsStrong.interpolate( real_c( -2.0 ) * ( COEFF_A20 + COEFF_A02 ), level );

   // Run tests
   WALBERLA_LOG_INFO_ON_ROOT( "Running tests with (P2ConstantMassOperator, P2ConstantLaplaceOperator)" );
   testSmoother<P2ConstantMassOperator, P2ConstantLaplaceOperator>( GAUSS_SEIDEL, rhsStrong, poly, storage, level );
   testSmoother<P2ConstantMassOperator, P2ConstantLaplaceOperator>( JACOBI, rhsStrong, poly, storage, level );

   WALBERLA_LOG_INFO_ON_ROOT( "Running tests with (P2ElementwiseMassOperator, P2ConstantLaplaceOperator)" );
   testSmoother<P2ElementwiseMassOperator, P2ConstantLaplaceOperator>( JACOBI, rhsStrong, poly, storage, level );
   // testSmoother<P2ElementwiseMassOperator, P2ConstantLaplaceOperator>( GAUSS_SEIDEL, rhsStrong, poly, storage, level );

   WALBERLA_LOG_INFO_ON_ROOT( "Running tests with (P2ElementwiseMassOperator, P2ElementwiseLaplaceOperator)" );
   testSmoother<P2ElementwiseMassOperator, P2ElementwiseLaplaceOperator>( JACOBI, rhsStrong, poly, storage, level );

   return EXIT_SUCCESS;
}
