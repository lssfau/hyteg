/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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
#include "core/logging/Logging.h"
#include "core/timing/TimingJSON.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

uint_t   level;
uint_t   numProc;
MeshInfo meshInfo    = MeshInfo::emptyMeshInfo();
bool     printTiming = false;
bool     writeVTK    = false;
bool printCGInfo = false;

template < typename LaplaceOperator, typename MassOperator >
int runBenchmark( const double tolerance, const uint_t maxIter )
{
   LIKWID_MARKER_REGISTER( "solve" );
   SetupPrimitiveStorage setupStorage( meshInfo, numProc );
   //auto storage = PrimitiveStorage( setupStorage );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );
   //auto storage = PrimitiveStorage::createFromGmshFile( meshFile );

   typename LaplaceOperator::srcType r( "r", storage, level, level );
   typename LaplaceOperator::srcType f( "f", storage, level, level );
   typename LaplaceOperator::srcType u( "u", storage, level, level );
   typename LaplaceOperator::srcType u_exact( "u_exact", storage, level, level );
   typename LaplaceOperator::srcType err( "err", storage, level, level );
   typename LaplaceOperator::srcType npoints_helper( "npoints_helper", storage, level, level );

   MassOperator    M( storage, level, level );
   LaplaceOperator L( storage, level, level );

   std::function< double( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& x ) {
      return ( 1.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };
   std::function< double( const hyteg::Point3D& ) > rhs = []( const hyteg::Point3D& x ) {
      return ( 3.0 / 2.0 ) * sin( 2 * x[0] ) * sinh( x[1] );
   };

   u.interpolate( exact, level, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, level );
   npoints_helper.interpolate( rhs, level );
   M.apply( npoints_helper, f, level, hyteg::All );

   auto solver = hyteg::CGSolver< LaplaceOperator >( storage, level, level, maxIter, tolerance );

   solver.setPrintInfo( printCGInfo );

   auto globalInfo = storage->getGlobalInfo();
   WALBERLA_LOG_INFO_ON_ROOT( globalInfo )
   printFunctionAllocationInfo( *storage, 1 );

   LIKWID_MARKER_START( "solve" );
   solver.solve( L, u, f, level );
   LIKWID_MARKER_STOP( "solve" );

   err.assign( {1.0, -1.0}, {u, u_exact}, level );

   const uint_t totalDoFs    = numberOfGlobalDoFs< typename LaplaceOperator::Operator::srcType::Tag >( *storage, level );
   const double discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / real_c( totalDoFs ) );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );

   if ( printTiming )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "creating timing tree" );
      walberla::WcTimingTree tt     = storage->getTimingTree()->getReduced().getCopyWithRemainder();
      nlohmann::json         ttjson = nlohmann::json( tt );
      std::ofstream          o( "P1CGBenchmarkOutput.json" );
      o << ttjson;
      o.close();
      WALBERLA_LOG_INFO_ON_ROOT( tt );
   }

   if ( writeVTK )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Writing VTK" );
      VTKOutput vtkOutput( "output", "cg_benchmark", storage );
      vtkOutput.add( u );
      vtkOutput.add( u_exact );
      vtkOutput.add( err );
      vtkOutput.add( f );
      vtkOutput.add( r );
      vtkOutput.write( level );
   }

   return 0;
}

int main( int argc, char* argv[] )
{
   LIKWID_MARKER_INIT;
   walberla::Environment env( argc, argv );
   //walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();
   LIKWID_MARKER_THREADINIT;

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./P1CGBenchmark.prm";
      WALBERLA_LOG_INFO_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   //Main Config
   level                            = mainConf.getParameter< uint_t >( "level" );
   const std::string discretization = mainConf.getParameter< std::string >( "discretization" );
   const std::string dimension      = mainConf.getParameter< std::string >( "dimension" );
   writeVTK                         = mainConf.getParameter< bool >( "writeVTK" );
   printTiming                      = mainConf.getParameter< bool >( "printTiming" );
   printCGInfo                      = mainConf.getParameter< bool >( "printCGInfo" );

   //2D Config
   const uint_t facesPerProcess = mainConf.getParameter< uint_t >( "facesPerProcess" );

   //3D Config
   const uint_t cubesX = mainConf.getParameter< uint_t >( "cubesX" );
   const uint_t cubesY = mainConf.getParameter< uint_t >( "cubesY" );
   const uint_t cubesZ = mainConf.getParameter< uint_t >( "cubesZ" );

   const double tolerance = 1e-15;
   const uint_t maxIter   = 1000000;

   numProc = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );
   if ( dimension == "2D" )
   {
      meshInfo = hyteg::MeshInfo::meshFaceChain( numProc * facesPerProcess );
   }
   else if ( dimension == "3D" )
   {
      meshInfo = hyteg::MeshInfo::meshSymmetricCuboid( Point3D( 0, 0, 0 ), Point3D( 1, 1, 1 ), cubesX, cubesY, cubesZ );
   }
   else
   {
      WALBERLA_ABORT( "Wrong dimension: " << dimension )
   }

   if ( discretization == "P1" )
   {
      runBenchmark< hyteg::P1ConstantLaplaceOperator, hyteg::P1ConstantMassOperator >( tolerance, maxIter );
   }
   else if ( discretization == "P2" )
   {
      runBenchmark< hyteg::P2ConstantLaplaceOperator, hyteg::P2ConstantMassOperator >( tolerance, maxIter );
   }
   else
   {
      WALBERLA_ABORT( "Wrong discretization: " << discretization )
   }
   LIKWID_MARKER_CLOSE;
}
