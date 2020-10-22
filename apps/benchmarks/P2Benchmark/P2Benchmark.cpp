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

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/misc/dummy.hpp"

using walberla::real_c;
using walberla::real_t;
using namespace hyteg;

/*
 * This benchmark meassures the time for several P2 functions on a macro face
 */
int main( int argc, char** argv )
{
   LIKWID_MARKER_INIT;

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   auto timingTree = std::make_shared< walberla::WcTimingTree >();
   walberla::WcTimer timer;

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared<walberla::config::Config>();
   if( env.config() == nullptr ) {
      auto defaultFile = "./P2Benchmark.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT("No Parameter file given loading default parameter file: " << defaultFile);
      cfg->readParameterFile( defaultFile );
   } else {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const uint_t level         = mainConf.getParameter< uint_t >( "level" );
   const std::string meshFile = mainConf.getParameter< std::string >( "mesh" );

   LIKWID_MARKER_THREADINIT;

   MeshInfo                            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage               setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 0, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   auto storageInfo = storage->getGlobalInfo();
   auto numVertexDoFs = numberOfGlobalDoFs< VertexDoFFunctionTag >( *storage, level );
   auto numEdgeDoFs   = numberOfGlobalDoFs< EdgeDoFFunctionTag >( *storage, level );
   auto numP2DoFsTotal = numberOfGlobalDoFs< P2FunctionTag >( *storage, level );

   WALBERLA_LOG_DEVEL_ON_ROOT( "" );
   WALBERLA_LOG_DEVEL_ON_ROOT( "=================================" );
   WALBERLA_LOG_DEVEL_ON_ROOT( "===== P2 Function Benchmark =====" );
   WALBERLA_LOG_DEVEL_ON_ROOT( "=================================" );
   WALBERLA_LOG_DEVEL_ON_ROOT( "" );

   WALBERLA_LOG_INFO_ON_ROOT( storageInfo );
   WALBERLA_LOG_INFO_ON_ROOT( "mesh:             " << meshFile );
   WALBERLA_LOG_INFO_ON_ROOT( "refinement level: " << level );
   WALBERLA_LOG_INFO_ON_ROOT( "# vertexdofs:     " << numVertexDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( "# edgedofs:       " << numEdgeDoFs );
   WALBERLA_LOG_INFO_ON_ROOT( "# total dofs:     " << numP2DoFsTotal );

   WALBERLA_LOG_INFO_ON_ROOT( "=== Starting measurements ===" );

   P2Function< real_t > src( "src", storage, level, level );
   P2Function< real_t > dst( "dst", storage, level, level );

   hyteg::P2ConstantLaplaceOperator M( storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > someFunction = [&]( const hyteg::Point3D& point )
   {
      return point[0] + point[1];
   };

   const uint_t numIterations = 3;

   LIKWID_MARKER_START( "interpolate constant" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      src.interpolate( 42.0, level );
   timer.end();
   LIKWID_MARKER_STOP( "interpolate constant" );
   WALBERLA_LOG_INFO_ON_ROOT( "interpolate constant: " << timer.last() );

   LIKWID_MARKER_START( "interpolate function" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      src.interpolate( someFunction, level );
   timer.end();
   LIKWID_MARKER_STOP( "interpolate function" );
   WALBERLA_LOG_INFO_ON_ROOT( "interpolate function: " << timer.last() );

   LIKWID_MARKER_START( "assign" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      dst.assign( {1.0}, {src}, level );
   timer.end();
   LIKWID_MARKER_STOP( "assign" );
   WALBERLA_LOG_INFO_ON_ROOT( "assign:               " << timer.last() );

   LIKWID_MARKER_START( "assign scaled" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      dst.assign( {1.23}, {src}, level );
   timer.end();
   LIKWID_MARKER_STOP( "assign scaled" );
   WALBERLA_LOG_INFO_ON_ROOT( "assign scaled:        " << timer.last() );

   LIKWID_MARKER_START( "apply" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      M.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply:                " << timer.last() );

   LIKWID_MARKER_START( "SOR" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      M.smooth_sor( src, dst, 1.1, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "SOR" );
   WALBERLA_LOG_INFO_ON_ROOT( "SOR:                  " << timer.last() );

   LIKWID_MARKER_START( "dot" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      src.dotGlobal( dst, level );
   timer.end();
   LIKWID_MARKER_STOP( "dot" );
   WALBERLA_LOG_INFO_ON_ROOT( "dot:                  " << timer.last() );

   LIKWID_MARKER_START( "dot self" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      dst.dotGlobal( dst, level );
   timer.end();
   LIKWID_MARKER_STOP( "dot self" );
   WALBERLA_LOG_INFO_ON_ROOT( "dot self:             " << timer.last() );

   LIKWID_MARKER_START( "sync all" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      communication::syncP2FunctionBetweenPrimitives( dst, level );
   timer.end();
   LIKWID_MARKER_STOP( "sync all" );
   WALBERLA_LOG_INFO_ON_ROOT( "sync all:             " << timer.last() );

   misc::dummy( &dst );
   misc::dummy( &src );

   auto timingTreeReducedWithRemainder = timingTree->getReduced().getCopyWithRemainder();
   WALBERLA_LOG_INFO_ON_ROOT( timingTreeReducedWithRemainder );

   nlohmann::json ttjson = nlohmann::json( timingTreeReducedWithRemainder );
   std::ofstream o("P2Benchmark.json");
   o << ttjson;
   o.close();

   LIKWID_MARKER_CLOSE;
}
