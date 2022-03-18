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

// Test is to ensure that we can solve a small 2D problem for
// VectorLaplace, Epsilon and FullViscous operators with CG

#include "core/Environment.h"
#include "core/logging/Logging.h"
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
#include "hyteg/solvers/CGSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

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

   typename opType::srcType u( "u", storage, level, level );
   typename opType::srcType f( "f", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > guessX = []( const hyteg::Point3D& x ) {
      return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };

   std::function< real_t( const hyteg::Point3D& ) > guessY = []( const hyteg::Point3D& x ) {
      return ( 1.0 / 2.0 ) * sin( 2 * x[1] ) * sinh( x[2] );
   };

   u.interpolate( {guessX, guessY}, level );
   u.interpolate( real_c( 0 ), level, DirichletBoundary );
   f.interpolate( real_c( 0 ), level );

   uint_t maxIter  = 150;
   real_t tol      = 1e-10;
   auto   cgSolver = hyteg::CGSolver< opType >( storage, level, level, maxIter, tol );
   cgSolver.setPrintInfo( verbose );

   cgSolver.solve( *oper, u, f, level );

   typename FunctionTrait< typename opType::srcType >::VectorComponentType tmp1( "tmp1", storage, level, level );
   typename FunctionTrait< typename opType::srcType >::VectorComponentType tmp2( "tmp2", storage, level, level );

   real_t error = velocityMaxMagnitude( u, tmp1, tmp2, level, All );
   WALBERLA_LOG_INFO_ON_ROOT( "   --> error = " << error );

   // output data for visualisation
   if ( verbose )
   {
      std::string fName = "VectorToVectorOperatorCGTest_for_" + opName;
      VTKOutput   vtkOutput( "../../output", fName, storage );
      vtkOutput.add( u );
      vtkOutput.add( f );
      vtkOutput.write( level );
   }

   // Check that CG converged
   WALBERLA_CHECK_LESS( error, tol );
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( "----------------------------\n"
                              << " Testing application of CG\n"
                              << "----------------------------" );

   runCheck< P1ConstantVectorLaplaceOperator >( "P1ConstantVectorLaplaceOperator" );
   runCheck< P2ConstantVectorLaplaceOperator >( "P2ConstantVectorLaplaceOperator" );

   runCheck< P1ConstantEpsilonOperator >( "P1ConstantEpsilonOperator" );
   runCheck< P2ConstantEpsilonOperator >( "P2ConstantEpsilonOperator" );

   runCheck< P1ElementwiseAffineEpsilonOperator, true >( "P1ElementwiseAffineEpsilonOperator" );
   runCheck< P2ElementwiseAffineEpsilonOperator, true >( "P2ElementwiseAffineEpsilonOperator" );

   runCheck< P1ElementwiseBlendingEpsilonOperator, true >( "P1ElementwiseBlendingEpsilonOperator" );
   runCheck< P2ElementwiseBlendingEpsilonOperator, true >( "P2ElementwiseBlendingEpsilonOperator" );

   runCheck< P2ConstantFullViscousOperator >( "P2ConstantFullViscousOperator" );
   runCheck< P2ElementwiseBlendingFullViscousOperator, true >( "P2ElementwiseBlendingFullViscousOperator" );
}
