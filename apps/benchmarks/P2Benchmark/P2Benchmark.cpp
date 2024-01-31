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
#include "hyteg/forms/form_hyteg_generated/p2/p2_mass_blending_q4.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/geometry/IdentityMap.hpp"
#include "hyteg/misc/dummy.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "constantStencilOperator/P2ConstantOperator.hpp"

using walberla::real_c;
using walberla::real_t;
using namespace hyteg;

/*
 * This benchmark meassures the time for several P2 functions on a macro face
 */

void runFunctionTests( std::shared_ptr< PrimitiveStorage >& storage, uint_t level,
                       const uint_t numIterations, walberla::WcTimer& timer )
{
   P2Function< real_t > src( "src", storage, level, level );
   P2Function< real_t > dst( "dst", storage, level, level );

   std::function< real_t( const hyteg::Point3D& ) > someFunction = [&]( const hyteg::Point3D& point ) {
     return point[0] + point[1];
   };

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
      communication::syncFunctionBetweenPrimitives( dst, level );
   timer.end();
   LIKWID_MARKER_STOP( "sync all" );
   WALBERLA_LOG_INFO_ON_ROOT( "sync all:                           " << timer.last() );


   misc::dummy( &dst );
   misc::dummy( &src );
}

void runLaplaceOperatorTests( SetupPrimitiveStorage& setupStorage, uint_t level,
                              const uint_t numIterations, walberla::WcTimer& timer,
                              std::shared_ptr< walberla::WcTimingTree >& timingTree )
{

   WALBERLA_LOG_INFO_ON_ROOT( "--- Timing Laplace.apply() ---" );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );
   bool running3D = storage->hasGlobalCells();

   WALBERLA_LOG_INFO_ON_ROOT( "Allocating functions ..." );

   P2Function< real_t > src( "src", storage, level, level );
   P2Function< real_t > dst( "dst", storage, level, level );

   WALBERLA_LOG_INFO_ON_ROOT( "Preparing operators ..." );

   WALBERLA_LOG_INFO_ON_ROOT( "- L_constant_stencil ..." );
   P2ConstantLaplaceOperator L_constant_stencil( storage, level, level );

   WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_otf_cc ..." );
   P2ElementwiseLaplaceOperator L_elementwise_otf_cc( storage, level, level );

   WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_otf_blending_id ..." );
   P2ElementwiseBlendingLaplaceOperator L_elementwise_otf_blending_id( storage, level, level, P2Form_laplace(), false );

   WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_stored_blending_id ..." );
   P2ElementwiseBlendingLaplaceOperator L_elementwise_stored_blending_id( storage, level, level, P2Form_laplace(), false );
   L_elementwise_stored_blending_id.computeAndStoreLocalElementMatrices();

   if ( running3D )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_otf_blending_shell ..." );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_otf_blending_annulus ..." );
   }
   P2ElementwiseBlendingLaplaceOperator L_elementwise_otf_blending_shell( storage, level, level, P2Form_laplace(), false );

   WALBERLA_LOG_INFO_ON_ROOT( "Done." );

   LIKWID_MARKER_START( "apply cc stencil" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_constant_stencil.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply cc stencil" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply cc stencil:              " << timer.last() );

   LIKWID_MARKER_START( "apply cc elem" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_elementwise_otf_cc.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply cc elem" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply cc elem:                 " << timer.last() );

   LIKWID_MARKER_START( "apply blending id elem" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_elementwise_otf_blending_id.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply blending id elem" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply blending id elem:        " << timer.last() );

   LIKWID_MARKER_START( "apply blending id elem stored" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_elementwise_stored_blending_id.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply blending id elem stored" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply blending id elem stored: " << timer.last() );

   if ( running3D )
   {
     IcosahedralShellMap::setMap( setupStorage );
     storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );
     LIKWID_MARKER_START( "apply blending shell elem" );
   }
   else
   {
     AnnulusMap::setMap( setupStorage );
     storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );
     LIKWID_MARKER_START( "apply blending annulus elem" );
   }

   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_elementwise_otf_blending_shell.apply( src, dst, level, All );
   timer.end();
   if ( running3D )
   {
     LIKWID_MARKER_STOP( "apply blending shell elem" );
     WALBERLA_LOG_INFO_ON_ROOT( "apply blending shell elem:     " << timer.last() );
   }
   else
   {
     LIKWID_MARKER_STOP( "apply blending annulus elem" );
     WALBERLA_LOG_INFO_ON_ROOT( "apply blending annulus elem:     " << timer.last() );
   }

}


void runMassOperatorTests( SetupPrimitiveStorage& setupStorage, uint_t level,
                              const uint_t numIterations, walberla::WcTimer& timer,
                              std::shared_ptr< walberla::WcTimingTree >& timingTree )
{
   WALBERLA_LOG_INFO_ON_ROOT( "--- Timing Mass.apply() ---" );

   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );
   bool running3D = storage->hasGlobalCells();

   WALBERLA_LOG_INFO_ON_ROOT( "Allocating functions ..." );

   P2Function< real_t > src( "src", storage, level, level );
   P2Function< real_t > dst( "dst", storage, level, level );

   WALBERLA_LOG_INFO_ON_ROOT( "Preparing operators ..." );

   WALBERLA_LOG_INFO_ON_ROOT( "- L_constant_stencil ..." );
   P2ConstantMassOperator L_constant_stencil( storage, level, level );

   WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_otf_cc ..." );
   P2ElementwiseMassOperator L_elementwise_otf_cc( storage, level, level );

   WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_otf_blending_id ..." );
   P2ElementwiseBlendingMassOperator L_elementwise_otf_blending_id( storage, level, level, forms::p2_mass_blending_q5(), false );

   WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_stored_blending_id ..." );
   P2ElementwiseBlendingMassOperator L_elementwise_stored_blending_id( storage, level, level, forms::p2_mass_blending_q5(), false );
   L_elementwise_stored_blending_id.computeAndStoreLocalElementMatrices();

   if ( running3D )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_otf_blending_shell ..." );
   }
   else
   {
      WALBERLA_LOG_INFO_ON_ROOT( "- L_elementwise_otf_blending_annulus ..." );
   }
   P2ElementwiseBlendingMassOperator L_elementwise_otf_blending_shell( storage, level, level, forms::p2_mass_blending_q5(), false );

   WALBERLA_LOG_INFO_ON_ROOT( "Done." );

   LIKWID_MARKER_START( "apply cc stencil" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_constant_stencil.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply cc stencil" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply cc stencil:              " << timer.last() );

   LIKWID_MARKER_START( "apply cc elem" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_elementwise_otf_cc.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply cc elem" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply cc elem:                 " << timer.last() );

   LIKWID_MARKER_START( "apply blending id elem" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_elementwise_otf_blending_id.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply blending id elem" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply blending id elem:        " << timer.last() );

   LIKWID_MARKER_START( "apply blending id elem stored" );
   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_elementwise_stored_blending_id.apply( src, dst, level, All );
   timer.end();
   LIKWID_MARKER_STOP( "apply blending id elem stored" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply blending id elem stored: " << timer.last() );

   if ( running3D )
   {
     IcosahedralShellMap::setMap( setupStorage );
     storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );
     LIKWID_MARKER_START( "apply blending shell elem" );
   }
   else
   {
     AnnulusMap::setMap( setupStorage );
     storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );
     LIKWID_MARKER_START( "apply blending annulus elem" );
   }

   timer.reset();
   for ( uint_t i = 0; i < numIterations; i++ )
      L_elementwise_otf_blending_shell.apply( src, dst, level, All );
   timer.end();
   if ( running3D )
   {
     LIKWID_MARKER_STOP( "apply blending shell elem" );
     WALBERLA_LOG_INFO_ON_ROOT( "apply blending shell elem:     " << timer.last() );
   }
   else
   {
     LIKWID_MARKER_STOP( "apply blending annulus elem" );
     WALBERLA_LOG_INFO_ON_ROOT( "apply blending annulus elem:     " << timer.last() );
   }
}


void performBenchmarkRuns( const walberla::Config::BlockHandle& conf )
{
   auto              timingTree = std::make_shared< walberla::WcTimingTree >();
   walberla::WcTimer timer;

   const uint_t      level    = conf.getParameter< uint_t >( "level" );
   const uint_t      numIterations = conf.getParameter< uint_t >( "numIterations" );
   const std::string meshFile = conf.getParameter< std::string >( "mesh" );

   MeshInfo              meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 0, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   auto storageInfo    = storage->getGlobalInfo();
   auto numVertexDoFs  = numberOfGlobalDoFs< VertexDoFFunctionTag >( *storage, level );
   auto numEdgeDoFs    = numberOfGlobalDoFs< EdgeDoFFunctionTag >( *storage, level );
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
   WALBERLA_LOG_INFO_ON_ROOT( "numIterations:    " << numIterations );

   WALBERLA_LOG_INFO_ON_ROOT( "=== Starting measurements ===" );

   runFunctionTests( storage, level, numIterations, timer );
   if ( conf.getParameter< bool >( "timeLaplace" ) ) runLaplaceOperatorTests( setupStorage, level, numIterations, timer, timingTree );
   if ( conf.getParameter< bool >( "timeMass" ) ) runMassOperatorTests( setupStorage, level, numIterations, timer, timingTree );

   auto timingTreeReducedWithRemainder = timingTree->getReduced().getCopyWithRemainder();
   if ( conf.getParameter< bool >( "printDetails" ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( timingTreeReducedWithRemainder );
   }
   nlohmann::json ttjson = nlohmann::json( timingTreeReducedWithRemainder );
   const std::string jsonFile = conf.getParameter< std::string >( "jsonFile" );
   std::ofstream  o( jsonFile );
   o << ttjson;
   o.close();
}


int main( int argc, char** argv )
{
   LIKWID_MARKER_INIT;

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared<walberla::config::Config>();
   if( env.config() == nullptr ) {
      auto defaultFile = "./P2Benchmark.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT("No Parameter file given loading default parameter file: " << defaultFile);
      cfg->readParameterFile( defaultFile );
   } else {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle conf2D = cfg->getBlock( "Parameters_for_2D" );
   const walberla::Config::BlockHandle conf3D = cfg->getBlock( "Parameters_for_3D" );

   LIKWID_MARKER_THREADINIT;

   if ( conf2D.getParameter< bool >( "run" ) )
   {
     performBenchmarkRuns( conf2D );
   }

   if ( conf3D.getParameter< bool >( "run" ) )
   {
     performBenchmarkRuns( conf3D );
   }

   LIKWID_MARKER_CLOSE;
}
