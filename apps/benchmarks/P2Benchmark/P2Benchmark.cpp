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

#include <memory>

#include "core/Environment.h"
#include "core/timing/TimingJSON.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/edgedofspace/EdgeDoFMacroFace.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/IdentityMap.hpp"
#include "hyteg/misc/dummy.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

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

   auto              timingTree = std::make_shared< walberla::WcTimingTree >();
   walberla::WcTimer timer;

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./P2Benchmark.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const uint_t      level    = mainConf.getParameter< uint_t >( "level" );
   const std::string meshFile = mainConf.getParameter< std::string >( "mesh" );

   LIKWID_MARKER_THREADINIT;

   MeshInfo meshInfo = MeshInfo::fromGmshFile( meshFile );
   auto     setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 0, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( *setupStorage, timingTree );

   auto storageInfo    = storage->getGlobalInfo();
   auto numVertexDoFs  = numberOfGlobalDoFs< VertexDoFFunctionTag >( *storage, level );
   auto numEdgeDoFs    = numberOfGlobalDoFs< EdgeDoFFunctionTag >( *storage, level );
   auto numP2DoFsTotal = numberOfGlobalDoFs< P2FunctionTag >( *storage, level );

   const uint_t numIterations = 3;

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
   WALBERLA_LOG_INFO_ON_ROOT( "# iterations:     " << numIterations );

   WALBERLA_LOG_INFO_ON_ROOT( "Allocating functions ..." );

   P2Function< real_t > src( "src", storage, level, level );
   P2Function< real_t > dst( "dst", storage, level, level );

   WALBERLA_LOG_INFO_ON_ROOT( "Preparing operators ..." );

   WALBERLA_LOG_INFO_ON_ROOT( "- L_constant_stencil ..." );
   P2ConstantLaplaceOperator            L_constant_stencil( storage, level, level );
   WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_otf_cc ..." );
   P2ElementwiseLaplaceOperator         L_elementwise_otf_cc( storage, level, level );
   WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_otf_blending_id ..." );
   P2ElementwiseBlendingLaplaceOperator L_elementwise_otf_blending_id( storage, level, level, P2Form_laplace(), false );
   WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_otf_blending_shell ..." );
   P2ElementwiseBlendingLaplaceOperator L_elementwise_otf_blending_shell( storage, level, level, P2Form_laplace(), false );
   WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_stored_blending_shell ..." );
   P2ElementwiseBlendingLaplaceOperator L_elementwise_stored_blending_shell( storage, level, level, P2Form_laplace(), false );

   WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_otf_blending_id_pimped ..." );
   P2ElementwiseBlendingLaplaceOperatorPimped3D L_elementwise_otf_blending_id_pimped( storage, level, level, P2Form_laplacePimped3D(), false );
   WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_otf_blending_shell_pimped ..." );
   P2ElementwiseBlendingLaplaceOperatorPimped3D L_elementwise_otf_blending_shell_pimped( storage, level, level, P2Form_laplacePimped3D(), false );

   WALBERLA_LOG_INFO_ON_ROOT( "Done." );

   std::function< real_t( const hyteg::Point3D& ) > someFunction = [&]( const hyteg::Point3D& point ) {
     return point[0] + point[1];
   };

   WALBERLA_LOG_INFO_ON_ROOT( "=== Starting measurements ===" );

   LIKWID_MARKER_START( "interpolate constant" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      src.interpolate( 42.0, level );
   timer.end();
   LIKWID_MARKER_STOP( "interpolate constant" );
   WALBERLA_LOG_INFO_ON_ROOT( "interpolate constant:               " << timer.last() );

   LIKWID_MARKER_START( "interpolate function" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      src.interpolate( someFunction, level );
   timer.end();
   LIKWID_MARKER_STOP( "interpolate function" );
   WALBERLA_LOG_INFO_ON_ROOT( "interpolate function:               " << timer.last() );

   LIKWID_MARKER_START( "assign" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      dst.assign( {1.0}, {src}, level );
   timer.end();
   LIKWID_MARKER_STOP( "assign" );
   WALBERLA_LOG_INFO_ON_ROOT( "assign:                             " << timer.last() );

   LIKWID_MARKER_START( "assign scaled" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      dst.assign( {1.23}, {src}, level );
   timer.end();
   LIKWID_MARKER_STOP( "assign scaled" );
   WALBERLA_LOG_INFO_ON_ROOT( "assign scaled:                      " << timer.last() );

   LIKWID_MARKER_START( "dot" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      src.dotGlobal( dst, level );
   timer.end();
   LIKWID_MARKER_STOP( "dot" );
   WALBERLA_LOG_INFO_ON_ROOT( "dot:                                " << timer.last() );

   LIKWID_MARKER_START( "dot self" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      dst.dotGlobal( dst, level );
   timer.end();
   LIKWID_MARKER_STOP( "dot self" );
   WALBERLA_LOG_INFO_ON_ROOT( "dot self:                           " << timer.last() );

   LIKWID_MARKER_START( "sync all" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      communication::syncP2FunctionBetweenPrimitives( dst, level );
   timer.end();
   LIKWID_MARKER_STOP( "sync all" );
   WALBERLA_LOG_INFO_ON_ROOT( "sync all:                           " << timer.last() );

   LIKWID_MARKER_START( "apply cc stencil" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_constant_stencil.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply cc stencil" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply cc stencil:                   " << timer.last() );

   LIKWID_MARKER_START( "apply cc elem" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_elementwise_otf_cc.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply cc elem" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply cc elem:                      " << timer.last() );

   LIKWID_MARKER_START( "apply blending id elem" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_elementwise_otf_blending_id.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply blending id elem" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply blending id elem:             " << timer.last() );

   LIKWID_MARKER_START( "apply blending id elem (pimped)" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_elementwise_otf_blending_id_pimped.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply blending id elem (pimped)" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply blending id elem (pimped):    " << timer.last() );

   IcosahedralShellMap::setMap( *setupStorage );
   storage = std::make_shared< PrimitiveStorage >( *setupStorage, timingTree );

   LIKWID_MARKER_START( "apply blending shell elem" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_elementwise_otf_blending_shell.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply blending shell elem" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply blending shell elem:          " << timer.last() );

   LIKWID_MARKER_START( "apply blending shell elem (pimped)" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_elementwise_otf_blending_shell_pimped.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply blending shell elem (pimped)" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply blending shell elem (pimped): " << timer.last() );

   L_elementwise_stored_blending_shell.computeAndStoreLocalElementMatrices();
   LIKWID_MARKER_START( "apply blending shell elem stored" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_elementwise_stored_blending_shell.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply blending shell elem stored" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply blending shell elem stored:   " << timer.last() );
   LIKWID_MARKER_START( "apply blending shell elem" );

   meshInfo = MeshInfo::fromGmshFile( meshFile );
   setupStorage =
       std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage->setMeshBoundaryFlagsOnBoundary( 0, 0, true );
   storage = std::make_shared< PrimitiveStorage >( *setupStorage, timingTree );

   LIKWID_MARKER_START( "SOR cc stencil" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_constant_stencil.smooth_sor( src, dst, 1.1, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "SOR cc stencil" );
   WALBERLA_LOG_INFO_ON_ROOT( "SOR cc stencil:                     " << timer.last() );

   misc::dummy( &dst );
   misc::dummy( &src );

   auto timingTreeReducedWithRemainder = timingTree->getReduced().getCopyWithRemainder();
   WALBERLA_LOG_INFO_ON_ROOT( timingTreeReducedWithRemainder );

   nlohmann::json ttjson = nlohmann::json( timingTreeReducedWithRemainder );
   std::ofstream  o( "P2Benchmark.json" );
   o << ttjson;
   o.close();

   LIKWID_MARKER_CLOSE;
}
