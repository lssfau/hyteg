/*
 * Copyright (c) 2020 Andreas Wagner
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
#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/operators/VectorLaplaceOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1EpsilonOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"
#include "hyteg/p2functionspace/P2FullViscousOperator.hpp"
#include "hyteg/p2functionspace/P2VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/SORSmoother.hpp"
#include "hyteg/solvers/SymmetricGaussSeidelSmoother.hpp"
#include "hyteg/solvers/SymmetricSORSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"

using walberla::real_t;
using walberla::math::pi;

using namespace hyteg;

template < typename SmootherType, typename OperatorType, typename FunctionType >
bool smootherThrowsException( SmootherType& smoother, OperatorType& op, FunctionType& src, FunctionType& dst, uint_t level )
{
   try
   {
      smoother.solve( op, dst, src, level );
   } catch ( const std::runtime_error& e )
   {
      // WALBERLA_LOG_INFO_ON_ROOT( e.what() );
      return true;
   }
   return false;
}

template < typename opType, bool forceComputeInvDiag = true, bool needsViscosity = false >
void runCheck( const std::array< bool, 6 > properties, std::string opName )
{
   WALBERLA_LOG_INFO_ON_ROOT( "-----------------------------------------------\n"
                              << "Testing Properties of " << opName << "\n-----------------------------------------------" );

   bool wJacSmoothable = properties[0];
   bool gsSmoothable   = properties[1];
   bool sorSmoothable  = properties[2];
   bool sgsSmoothable  = properties[3];
   bool ssorSmoothable = properties[4];
   bool chebSmoothable = properties[5];

   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( {-1, -1} ), Point2D( {1., 1.} ), MeshInfo::CRISSCROSS, 2, 2 );
   // MeshInfo meshInfo = MeshInfo::meshCuboid( Point3D( {-1, -1, 0} ), Point3D( {1, 1, 2} ), 2, 2, 2 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage  = std::make_shared< PrimitiveStorage >( setupStorage );
   const uint_t                        minLevel = 2;
   const uint_t                        maxLevel = 3;

   using fType = typename opType::srcType;

   fType src( "srcFunc", storage, minLevel, maxLevel );
   fType dst( "dstFunc", storage, minLevel, maxLevel );

   std::shared_ptr< opType > oper;

   if constexpr ( needsViscosity )
   {
      std::function< real_t( const Point3D& ) > mu = []( const Point3D& x ) { return real_c( 3 ) * x[0] + x[1]; };

      oper = std::make_shared< opType >( storage, minLevel, maxLevel, mu );
   }
   else
   {
      oper = std::make_shared< opType >( storage, minLevel, maxLevel );
   }

   WeightedJacobiSmoother< opType >       wJacSmoother( storage, minLevel, maxLevel, real_c( 0.5 ) );
   GaussSeidelSmoother< opType >          gsSmoother;
   SymmetricGaussSeidelSmoother< opType > sgsSmoother;
   SORSmoother< opType >                  sorSmoother( real_c( 1.5 ) );
   SymmetricSORSmoother< opType >         ssorSmoother( real_c( 1.5 ) );
   ChebyshevSmoother< opType >            chebSmoother( storage, minLevel, maxLevel );

   // Chebyshev needs additional preparations
   if constexpr ( forceComputeInvDiag )
   {
      // WALBERLA_LOG_INFO_ON_ROOT( "Executing <operator>::computeInverseDiagonalOperatorValues" );
      oper->computeInverseDiagonalOperatorValues();
   }
   chebSmoother.setupCoefficients( 1, 1 );

   bool hasProperty = !smootherThrowsException( wJacSmoother, *oper, src, dst, minLevel );
   WALBERLA_CHECK( hasProperty == wJacSmoothable );
   WALBERLA_LOG_INFO_ON_ROOT( "* weighted Jacobi smoothable .......... " << ( hasProperty ? "yes" : "no" ) );

   hasProperty = !smootherThrowsException( gsSmoother, *oper, src, dst, minLevel );
   WALBERLA_CHECK( hasProperty == gsSmoothable );
   WALBERLA_LOG_INFO_ON_ROOT( "* Gauss-Seidel smoothable ............. " << ( hasProperty ? "yes" : "no" ) );

   hasProperty = !smootherThrowsException( sorSmoother, *oper, src, dst, minLevel );
   WALBERLA_CHECK( hasProperty == sorSmoothable );
   WALBERLA_LOG_INFO_ON_ROOT( "* SOR smoothable ...................... " << ( hasProperty ? "yes" : "no" ) );

   hasProperty = !smootherThrowsException( sgsSmoother, *oper, src, dst, minLevel );
   WALBERLA_CHECK( hasProperty == sgsSmoothable );
   WALBERLA_LOG_INFO_ON_ROOT( "* symmetric Gauss-Seidel smoothable ... " << ( hasProperty ? "yes" : "no" ) );

   hasProperty = !smootherThrowsException( ssorSmoother, *oper, src, dst, minLevel );
   WALBERLA_CHECK( hasProperty == ssorSmoothable );
   WALBERLA_LOG_INFO_ON_ROOT( "* symmetric SOR smoothable ............ " << ( hasProperty ? "yes" : "no" ) );

   hasProperty = !smootherThrowsException( chebSmoother, *oper, src, dst, minLevel );
   WALBERLA_CHECK( hasProperty == chebSmoothable );
   WALBERLA_LOG_INFO_ON_ROOT( "* Chebyshev ........................... " << ( hasProperty ? "yes" : "no" ) );
}

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   // =================
   //  ScalarOperators
   // =================
   runCheck< P1ConstantLaplaceOperator >( {true, true, true, true, true, true}, "P1ConstantLaplaceOperator" );
   runCheck< P2ConstantLaplaceOperator, false >( {true, true, true, false, false, false}, "P2ConstantLaplaceOperator (in 2D)" );

   runCheck< P1BlendingLaplaceOperator >( {true, true, true, false, false, true}, "P1BlendingLaplaceOperator" );

   runCheck< P1ElementwiseLaplaceOperator, false >( {true, false, false, false, false, true}, "P1ElementwiseLaplaceOperator" );
   runCheck< P2ElementwiseLaplaceOperator, false >( {true, false, false, false, false, true}, "P2ElementwiseLaplaceOperator" );

   // =========================
   //  VectorToVectorOperators
   // =========================
   runCheck< P1ConstantVectorLaplaceOperator >( {true, true, true, true, true, true}, "P1ConstantVectorLaplaceOperator" );

   runCheck< P2ConstantVectorLaplaceOperator, false >( {true, true, true, false, false, false},
                                                       "P2ConstantVectorLaplaceOperator (in 2D)" );

   runCheck< P1ElementwiseAffineEpsilonOperator, true, true >( {false, false, false, false, false, true},
                                                               "P1ElementwiseAffineEpsilonOperator" );

   runCheck< P1ElementwiseBlendingEpsilonOperator, true, true >( {false, false, false, false, false, true},
                                                                 "P1ElementwiseBlendingEpsilonOperator" );

   runCheck< P2ElementwiseBlendingFullViscousOperator, true, true >( {false, false, false, false, false, true},
                                                                     "P2ElementwiseBlendingFullViscousOperator" );
}
