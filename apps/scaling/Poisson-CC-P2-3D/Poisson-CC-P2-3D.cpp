/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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
#include "core/timing/TimingJSON.h"
#include "core/timing/TimingTree.h"

#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticRestriction.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/petsc/PETScLUSolver.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"
#include "hyteg/solvers/GeometricMultigridSolver.hpp"
#include "hyteg/solvers/SORSmoother.hpp"

using walberla::real_c;
using walberla::real_t;

using walberla::uint_c;
using walberla::uint_t;

namespace hyteg {

std::shared_ptr< PrimitiveStorage >
    createDomain( std::shared_ptr< walberla::WcTimingTree > timingTree, uint_t nx, uint_t ny, uint_t nz )
{
   Point3D lowerLeftFront( 0, 0, 0 );
   Point3D upperRightBack( 1, 1, 1 );
   auto    meshInfo = MeshInfo::meshCuboid( lowerLeftFront, upperRightBack, nx, ny, nz );

   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );

   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   return storage;
}

class ErrorAnalysis
{
 public:
   ErrorAnalysis( const std::shared_ptr< PrimitiveStorage >& storage )
   : storage_( storage )
   , active_( false )
   , r( "r", storage )
   , uExact( "uExact", storage )
   , error( "error", storage )
   {}

   ErrorAnalysis( const std::shared_ptr< PrimitiveStorage >&     storage,
                  uint_t                                         minLevel,
                  uint_t                                         maxLevel,
                  std::function< real_t( const hyteg::Point3D& ) > analyticalSolution )
   : storage_( storage )
   , active_( true )
   , r( "r", storage, minLevel, maxLevel )
   , uExact( "uExact", storage, minLevel, maxLevel )
   , error( "error", storage, minLevel, maxLevel )
   {
      for ( uint_t level = minLevel; level <= maxLevel; level++ )
      {
         uExact.interpolate( analyticalSolution, level );
      }
   }

   void update( const P2ConstantLaplaceOperator& A, const P2Function< real_t >& u, const P2Function< real_t >& f, uint_t level )
   {
      if ( !active_ )
         return;

      discreteL2ErrorLast_    = discreteL2Error_;
      discreteL2ResidualLast_ = discreteL2Residual_;

      const real_t numDoFs = real_c( numberOfGlobalDoFs< P2FunctionTag >( *storage_, level ) );

      A.apply( u, error, level, hyteg::Inner );
      r.assign( {1.0, -1.0}, {f, error}, level, hyteg::Inner );
      error.assign( {1.0, -1.0}, {u, uExact}, level );

      discreteL2Error_    = std::sqrt( error.dotGlobal( error, level, DoFType::All ) / numDoFs );
      discreteL2Residual_ = std::sqrt( r.dotGlobal( r, level, DoFType::Inner ) / numDoFs );

      discreteL2ErrorRate_    = discreteL2Error_ / discreteL2ErrorLast_;
      discreteL2ResidualRate_ = discreteL2Residual_ / discreteL2ResidualLast_;
   }

   void printHeader()
   {
      WALBERLA_LOG_INFO_ON_ROOT(
          " After cycle... ||     L2 error | L2 error reduction ||  L2 residual | L2 residual reduction |" );
      WALBERLA_LOG_INFO_ON_ROOT(
          " ---------------++--------------+--------------------++--------------+-----------------------|" );
   }

   void print( uint_t numPreviousCycles )
   {
      if ( numPreviousCycles == 0 )
      {
         WALBERLA_LOG_INFO_ON_ROOT( "        initial || " << std::scientific << discreteL2Error_ << " |                --- || "
                                                          << discreteL2Residual_ << " |                   --- |" );
      }
      else
      {
         WALBERLA_LOG_INFO_ON_ROOT( " " << std::setw( 14 ) << numPreviousCycles << " || " << std::scientific << discreteL2Error_
                                        << " |       " << discreteL2ErrorRate_ << " || " << discreteL2Residual_ << " |          "
                                        << discreteL2ResidualRate_ << " |" );
      }
   }

   real_t discreteL2Error_;
   real_t discreteL2ErrorLast_;
   real_t discreteL2ErrorRate_;
   real_t discreteL2Residual_;
   real_t discreteL2ResidualLast_;
   real_t discreteL2ResidualRate_;

 private:
   std::shared_ptr< PrimitiveStorage > storage_;
   bool                                active_;

   P2Function< real_t > r;
   P2Function< real_t > uExact;
   P2Function< real_t > error;
};

void run( uint_t                                    nx,
          uint_t                                    ny,
          uint_t                                    nz,
          uint_t                                    minLevel,
          uint_t                                    maxLevel,
          real_t                                    sorOmega,
          uint_t                                    preSmoothing,
          uint_t                                    postSmoothing,
          uint_t                                    numCycles,
          bool                                      calcError,
          std::shared_ptr< walberla::WcTimingTree > timingTree,
          bool                                      writeVTK )
{
   //////////////////
   // Domain setup //
   //////////////////

   timingTree->start( "01_Domain_Setup" );

   auto storage = createDomain( timingTree, nx, ny, nz );

   auto storageInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( storageInfo );

   const auto numberOfGlobalP2DoFs = numberOfGlobalDoFs< P2FunctionTag >( *storage, maxLevel );
   WALBERLA_LOG_INFO_ON_ROOT( "Number of global P2 DoFs (incl. boundary, level " << maxLevel << "): " << numberOfGlobalP2DoFs )

   if ( writeVTK )
   {
      writeDomainPartitioningVTK( storage, "output", "domain_partitioning" );
   }

   timingTree->stop( "01_Domain_Setup" );

   ///////////////////
   // Problem setup //
   ///////////////////

   timingTree->start( "02_Problem_Setup" );

   std::function< real_t( const hyteg::Point3D& ) > analyticalSolution = []( const hyteg::Point3D& x ) {
      return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };

   std::function< real_t( const hyteg::Point3D& ) > analyticalRHS = []( const hyteg::Point3D& x ) {
      return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };

   timingTree->start( "Function allocation" );

   P2Function< real_t > u( "u", storage, minLevel, maxLevel );
   P2Function< real_t > f( "f", storage, minLevel, maxLevel );

   timingTree->stop( "Function allocation" );

   timingTree->start( "Operator assembly" );

   P2ConstantMassOperator    M( storage, maxLevel, maxLevel );
   P2ConstantLaplaceOperator A( storage, minLevel, maxLevel );

   timingTree->stop( "Operator assembly" );

   timingTree->start( "Problem discretization" );

   u.interpolate( analyticalRHS, maxLevel, All );
   M.apply( u, f, maxLevel, Inner );

   u.interpolate( real_c( 0 ), maxLevel, All );
   u.interpolate( analyticalSolution, maxLevel, DirichletBoundary );

   timingTree->stop( "Problem discretization" );

   timingTree->stop( "02_Problem_Setup" );

   ////////////
   // Solver //
   ////////////

   timingTree->start( "03_Solver" );

   timingTree->start( "Solver setup" );

   auto smoother         = std::make_shared< SORSmoother< P2ConstantLaplaceOperator > >( sorOmega );
   auto restriction      = std::make_shared< P2toP2QuadraticRestriction >();
   auto prolongation     = std::make_shared< P2toP2QuadraticProlongation >();
   auto coarseGridSolver = std::make_shared< PETScLUSolver< P2ConstantLaplaceOperator > >( storage, minLevel );
   auto gmgSolver        = std::make_shared< GeometricMultigridSolver< P2ConstantLaplaceOperator > >(
       storage, smoother, coarseGridSolver, restriction, prolongation, minLevel, maxLevel, preSmoothing, postSmoothing );

   auto errorAnalysis = std::make_shared< ErrorAnalysis >( storage );
   if ( calcError )
   {
      errorAnalysis = std::make_shared< ErrorAnalysis >( storage, minLevel, maxLevel, analyticalSolution );
      errorAnalysis->update( A, u, f, maxLevel );
      errorAnalysis->printHeader();
      errorAnalysis->print( 0 );
   }

   timingTree->stop( "Solver setup" );

   timingTree->start( "Solver run" );

   for ( uint_t cycle = 0; cycle < numCycles; cycle++ )
   {
      gmgSolver->solve( A, u, f, maxLevel );
      errorAnalysis->update( A, u, f, maxLevel );
      errorAnalysis->print( cycle + 1 );
   }

   timingTree->stop( "Solver run" );

   timingTree->stop( "03_Solver" );
}

} // namespace hyteg

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::PETScManager petscManager( &argc, &argv );

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./parameters.prm";
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const uint_t minLevel      = mainConf.getParameter< uint_t >( "minLevel" );
   const uint_t maxLevel      = mainConf.getParameter< uint_t >( "maxlevel" );
   const real_t sorOmega      = mainConf.getParameter< real_t >( "sorOmega" );
   const uint_t preSmoothing  = mainConf.getParameter< uint_t >( "preSmoothing" );
   const uint_t postSmoothing = mainConf.getParameter< uint_t >( "postSmoothing" );
   const uint_t numCycles     = mainConf.getParameter< uint_t >( "numCycles" );
   const bool   calcError     = mainConf.getParameter< bool >( "calcError" );
   const bool   writeVTK      = mainConf.getParameter< bool >( "writeVTK" );
   const uint_t nx            = mainConf.getParameter< uint_t >( "nx" );
   const uint_t ny            = mainConf.getParameter< uint_t >( "ny" );
   const uint_t nz            = mainConf.getParameter< uint_t >( "nz" );

   const bool printTiming = mainConf.getParameter< bool >( "printTiming" );
   const bool jsonTiming  = mainConf.getParameter< bool >( "jsonTiming" );

   auto timingTree = std::make_shared< walberla::WcTimingTree >();

   hyteg::run( nx, ny, nz, minLevel, maxLevel, sorOmega, preSmoothing, postSmoothing, numCycles, calcError, timingTree, writeVTK );

   auto ttReduced = timingTree->getReduced().getCopyWithRemainder();
   if ( printTiming )
   {
      WALBERLA_LOG_INFO_ON_ROOT( ttReduced );
   }
   if ( jsonTiming )
   {
      nlohmann::json ttJson;
      walberla::timing::to_json( ttJson, ttReduced );
      std::ofstream jsonOutput;
      jsonOutput.open( "timingtree.json" );
      jsonOutput << ttJson.dump( 4 );
      jsonOutput.close();
      WALBERLA_LOG_INFO_ON_ROOT( "Wrote timing tree to JSON file." )
   }

   return 0;
}