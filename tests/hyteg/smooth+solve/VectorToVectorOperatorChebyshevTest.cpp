/*
 * Copyright (c) 2022 Marcus Mohr.
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

// Test is to ensure that we can perform Chebyshev smoothing for a small
// 2D problem and different VectorLaplace, Epsilon and FullViscous
// operators.
// Note: This is only an application not a convergence test.

#include "core/Environment.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/numerictools/CFDHelpers.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/p1functionspace/P1EpsilonOperator.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/p2functionspace/P2EpsilonOperator.hpp"
#include "hyteg/p2functionspace/P2FullViscousOperator.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

using namespace hyteg;

template < typename opType, bool needsViscosity = false >
void runCheck( std::string opName, bool verbose = false )
{
   WALBERLA_LOG_INFO_ON_ROOT( " * " << opName );

   const uint_t level = 3;
   // const std::string meshFile = "../../data/meshes/penta_5el.msh";
   const std::string meshFile = "../../data/meshes/tri_2el.msh";

   auto meshInfo = MeshInfo::fromGmshFile( meshFile );
   auto setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   std::shared_ptr< opType > oper;
   if constexpr ( needsViscosity )
   {
      std::function< real_t( const Point3D& ) > mu = []( const Point3D& x ) { return real_c( 3 ) + std::abs( x[0] + x[1] ); };

      oper = std::make_shared< opType >( storage, level, level, mu );
   }
   else
   {
      oper = std::make_shared< opType >( storage, level, level );
   }

   // we need the diagonal values
   if ( auto* op = dynamic_cast< OperatorWithInverseDiagonal< typename opType::srcType >* >( oper.get() ) )
   {
      op->computeInverseDiagonalOperatorValues();
   }
   else
   {
      WALBERLA_ABORT( "Houston we have a problem here!" );
   }

   typename opType::srcType u0( "u0", storage, level, level );
   typename opType::srcType u( "u", storage, level, level );
   typename opType::srcType f( "f", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > guessX = []( const hyteg::Point3D& x ) {
      return std::sin( real_c( 4 ) * pi * x[0] ) * std::cos( real_c( 5 ) * pi * x[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > guessY = []( const hyteg::Point3D& x ) {
      WALBERLA_UNUSED( x );
      return real_c( 0 );
      // return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };

   u.interpolate( guessX, level );
   typename opType::srcType tmpVec( "tmpVec", storage, level, level );
   real_t                   spectralRadius = chebyshev::estimateRadius( *oper, level, 100, storage, u, tmpVec );
   if ( verbose )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "   --> spectral estimate ........... " << spectralRadius );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "   --> spectrum estimation ......... check" );
   }

   auto chebyshev = hyteg::ChebyshevSmoother< opType >( storage, level, level );
   WALBERLA_LOG_INFO_ON_ROOT( "   --> smoother generation ......... check" );

   chebyshev.setupCoefficients( 2, spectralRadius );
   WALBERLA_LOG_INFO_ON_ROOT( "   --> setting coefficients ........ check" );

   f.interpolate( real_c( 0 ), level );
   u.interpolate( {guessX, guessY}, level );
   u.interpolate( real_c( 0 ), level, DirichletBoundary );
   u0.assign( {1.0}, {u}, level );

   chebyshev.solve( *oper, u, f, level );
   chebyshev.solve( *oper, u, f, level );
   WALBERLA_LOG_INFO_ON_ROOT( "   --> performing iterations ....... check" );

   if ( verbose )
   {
      real_t error = velocityMaxMagnitude( u, tmpVec[0], tmpVec[1], level, All );
      WALBERLA_LOG_INFO_ON_ROOT( "   --> error = " << error );
   }

   // output data for visualisation
   if ( verbose )
   {
      std::string fName = "VectorToVectorOperatorChebshevTest_for_" + opName;
      VTKOutput   vtkOutput( "../../output", fName, storage );
      vtkOutput.add( u );
      vtkOutput.add( u0 );
      vtkOutput.add( f );
      vtkOutput.write( level );
   }
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------\n"
                              << " Testing Chebyshev Smoothing\n"
                              << "-----------------------------" );

   runCheck< P1ConstantVectorLaplaceOperator >( "P1ConstantVectorLaplaceOperator" );

   runCheck< P1ConstantEpsilonOperator >( "P1ConstantEpsilonOperator" );

   runCheck< P1ElementwiseAffineEpsilonOperator, true >( "P1ElementwiseAffineEpsilonOperator" );
   runCheck< P2ElementwiseAffineEpsilonOperator, true >( "P2ElementwiseAffineEpsilonOperator" );

   runCheck< P1ElementwiseBlendingEpsilonOperator, true >( "P1ElementwiseBlendingEpsilonOperator" );
   runCheck< P2ElementwiseBlendingEpsilonOperator, true >( "P2ElementwiseBlendingEpsilonOperator" );

   runCheck< P2ElementwiseBlendingFullViscousOperator, true >( "P2ElementwiseBlendingFullViscousOperator" );

   // P2Constant operators are currently not Chebyshev smoothable, as the cannot
   // explicitely return their inverse diagonal
   // runCheck< P2ConstantVectorLaplaceOperator >( "P2ConstantVectorLaplaceOperator" );
   // runCheck< P2ConstantEpsilonOperator >( "P2ConstantEpsilonOperator" );
   // runCheck< P2ConstantFullViscousOperator >( "P2ConstantFullViscousOperator" );
}
