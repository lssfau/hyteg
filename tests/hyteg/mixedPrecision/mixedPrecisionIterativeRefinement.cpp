/*
* Copyright (c) 2017-2022 Dominik Thoennes.
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
#include "core/extern/json.hpp"
#include "core/logging/Logging.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/controlflow/SolverLoop.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

typedef P1ConstantOperator< P1FenicsForm< p1_diffusion_cell_integral_0_otherwise, p1_tet_diffusion_cell_integral_0_otherwise >,
                            false,
                            false,
                            false,
                            float >
    P1ConstantLaplaceOperatorSP;

double computeError( uint_t                      level,
                     P1Function< double >&       err,
                     P1ConstantMassOperator&     M,
                     P1Function< double >&       Merr,
                     const P1Function< double >& u,
                     P1Function< double >&       u_exact )
{
   err.assign( { 1.0, -1.0 }, { u, u_exact }, level, Inner );
   M.apply( err, Merr, level, Inner );
   return std::sqrt( err.dotGlobal( Merr, level, Inner ) );
}

void solvePoisson( uint_t minLevel, uint_t maxLevel )
{
   auto meshInfo = MeshInfo::meshUnitSquare( 0 );

   std::function< double( const hyteg::Point3D& ) > function_exact = []( const hyteg::Point3D& x ) {
      return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };

   std::function< double( const hyteg::Point3D& ) > function_rhs = []( const hyteg::Point3D& x ) {
      return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };

   auto setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( *setupStorage );

   P1Function< double > u( "u", storage, minLevel, maxLevel );
   P1Function< double > u_exact( "u", storage, minLevel, maxLevel );
   P1Function< double > tmp( "tmp", storage, minLevel, maxLevel );
   P1Function< double > f( "f", storage, minLevel, maxLevel );
   P1Function< double > r( "r", storage, minLevel, maxLevel );

   P1Function< float > uSP( "uSP", storage, minLevel, maxLevel );
   P1Function< float > rSP( "rSP", storage, minLevel, maxLevel );
   P1Function< float > tmpSP( "tmp", storage, minLevel, maxLevel );

   P1Function< double > err( "err", storage, minLevel, maxLevel );
   P1Function< double > Merr( "Merr", storage, minLevel, maxLevel );

   walberla::timing::Timer< walberla::timing::WcPolicy > timer;

   P1ConstantLaplaceOperator   A( storage, minLevel, maxLevel );
   P1ConstantLaplaceOperatorSP ASP( storage, minLevel, maxLevel );
   P1ConstantMassOperator      M( storage, minLevel, maxLevel );

   std::function< double( const Point3D& ) > random = []( const Point3D& x ) { return walberla::math::realRandom(); };

   for ( uint_t l = minLevel; l <= maxLevel; l++ )
   {
      u_exact.interpolate( function_exact, l, All );
      tmp.interpolate( function_rhs, l, All );
      u.interpolate( function_exact, l, DirichletBoundary );
      u.interpolate( random, l, Inner );
      M.apply( tmp, f, l, Inner );
      A.apply( u, r, l, Inner, Replace );
   }

   auto coarseGridSolver     = std::make_shared< CGSolver< P1ConstantLaplaceOperator > >( storage, minLevel, maxLevel );
   auto smoother             = std::make_shared< GaussSeidelSmoother< P1ConstantLaplaceOperator > >();
   auto restrictionOperator  = std::make_shared< P1toP1LinearRestriction< double > >();
   auto prolongationOperator = std::make_shared< P1toP1LinearProlongation< double > >();
   auto gmgSolver            = std::make_shared< GeometricMultigridSolver< P1ConstantLaplaceOperator > >(
       storage, smoother, coarseGridSolver, restrictionOperator, prolongationOperator, minLevel, maxLevel );

   auto coarseGridSolverSP     = std::make_shared< CGSolver< P1ConstantLaplaceOperatorSP > >( storage, minLevel, maxLevel );
   auto smootherSP             = std::make_shared< GaussSeidelSmoother< P1ConstantLaplaceOperatorSP > >();
   auto restrictionOperatorSP  = std::make_shared< P1toP1LinearRestriction< float > >();
   auto prolongationOperatorSP = std::make_shared< P1toP1LinearProlongation< float > >();
   auto gmgSolverSP            = std::make_shared< GeometricMultigridSolver< P1ConstantLaplaceOperatorSP > >(
       storage, smootherSP, coarseGridSolverSP, restrictionOperatorSP, prolongationOperatorSP, minLevel, maxLevel );

   double oldError, uError, previousLevelError = 0;

   std::map< std::string, std::function< void( uint_t ) > > differentSolvers;
   std::map< std::string, double >                          timings;
   std::map< std::string, double >                          errors;

   differentSolvers["Plain GMG"] = [&]( uint_t solveOnLevel ) { gmgSolver->solve( A, u, f, solveOnLevel ); };
   //
   differentSolvers["Iterative Refinement DP"] = [&]( uint_t solveOnLevel ) {
      A.apply( u, tmp, solveOnLevel, All, Replace );
      r.assign( { 1.0, -1.0 }, { f, tmp }, solveOnLevel, All );
      tmp.setToZero( solveOnLevel );
      gmgSolver->solve( A, tmp, r, solveOnLevel );
      u.assign( { 1.0, 1.0 }, { u, tmp }, solveOnLevel, Inner );
   };

   differentSolvers["Iterative Refinement SP"] = [&]( uint_t solveOnLevel ) {
      A.apply( u, tmp, solveOnLevel, All, Replace );
      r.assign( { 1.0, -1.0 }, { f, tmp }, solveOnLevel, All );
      rSP.copyFrom( r, solveOnLevel );
      tmpSP.setToZero( solveOnLevel );
      timer.start();
      gmgSolverSP->solve( ASP, tmpSP, rSP, solveOnLevel );
      timer.end();
      tmp.copyFrom( tmpSP, solveOnLevel );
      u.assign( { 1.0, 1.0 }, { u, tmp }, solveOnLevel, Inner );
   };

   for ( const auto& solver : differentSolvers )
   {
      WALBERLA_LOG_INFO( solver.first )
      for ( uint_t l = minLevel; l <= maxLevel; l++ )
      {
         u.interpolate( random, l, Inner );
      }
      timer.reset();
      for ( uint_t levelToSolve = minLevel; levelToSolve <= maxLevel; ++levelToSolve )
      {
         tmp.interpolate( function_rhs, levelToSolve, All );
         M.apply( tmp, f, levelToSolve, Inner );
         u.interpolate( function_exact, levelToSolve, DirichletBoundary );

         oldError = computeError( levelToSolve, err, M, Merr, u, u_exact );
         for ( ;; )
         {
            timer.start();
            solver.second( levelToSolve );
            timer.end();
            uError = computeError( levelToSolve, err, M, Merr, u, u_exact );
            if ( uError / oldError > 0.9 )
            {
               break;
            }
            oldError = uError;
         }

         WALBERLA_LOG_INFO( "Level " << levelToSolve << " Error: " << uError
                                     << " Error reduction: " << previousLevelError / uError )
         if ( levelToSolve > 2 )
         {
            WALBERLA_CHECK_GREATER( previousLevelError / uError, 3.4 )
         }
         previousLevelError = uError;

         if ( levelToSolve == maxLevel )
         {
            errors[solver.first]  = uError;
            timings[solver.first] = timer.total();
            WALBERLA_LOG_INFO( "=============================================================" )
         }
      }
   }
   for ( const auto& solver : differentSolvers )
   {
      WALBERLA_LOG_INFO( walberla::format(
          "%23s | L2Error: %9.3e | Time: %9.3e", solver.first.c_str(), errors[solver.first], timings[solver.first] ) );
   }
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( errors["Plain GMG"], errors["Iterative Refinement DP"], 1e-9 )
   WALBERLA_CHECK_FLOAT_EQUAL_EPSILON( errors["Plain GMG"], errors["Iterative Refinement SP"], 1e-9 )
}
} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::solvePoisson( 2, 7 );

   return 0;
}