#include <iostream>
#include <vector>

#include "core/Environment.h"
#include "core/timing/Timer.h"

#include "tinyhhg_core/LikwidWrapper.hpp"
#include "tinyhhg_core/edgedofspace/generatedKernels/all.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/misc/dummy.hpp"
#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/generatedKernels/all.hpp"
#include "tinyhhg_core/p1functionspace/generatedKernels/all.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/PrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"

using walberla::real_t;
using walberla::uint_t;

namespace hhg {

enum KernelType
{
   APPLY_V_TO_V_REPLACE,
   APPLY_V_TO_V_ADD,

   APPLY_V_TO_E_REPLACE,
   APPLY_V_TO_E_ADD,

   APPLY_E_TO_V_REPLACE,
   APPLY_E_TO_V_ADD,

   APPLY_E_TO_E_REPLACE,
   APPLY_E_TO_E_ADD,

   SOR_P1_V,
   SOR_P2_V,
   SOR_P2_E
};

const std::map< std::string, KernelType > strToKernelType = {
    {"APPLY_V_TO_V_REPLACE", APPLY_V_TO_V_REPLACE},
    {"APPLY_V_TO_V_ADD", APPLY_V_TO_V_ADD},

    {"APPLY_V_TO_E_REPLACE", APPLY_V_TO_E_REPLACE},
    {"APPLY_V_TO_E_ADD", APPLY_V_TO_E_ADD},

    {"APPLY_E_TO_V_REPLACE", APPLY_E_TO_V_REPLACE},
    {"APPLY_E_TO_V_ADD", APPLY_E_TO_V_ADD},

    {"APPLY_E_TO_E_REPLACE", APPLY_E_TO_E_REPLACE},
    {"APPLY_E_TO_E_ADD", APPLY_E_TO_E_ADD},

    {"SOR_P1_V", SOR_P1_V},
    {"SOR_P2_V", SOR_P2_V},
    {"SOR_P2_E", SOR_P2_E},
};

void runBenchmark( uint_t level, real_t iterationMinTime, KernelType kernelType )
{
   walberla::WcTimer timer;

   const uint_t chunkIterations = 100;

   const MeshInfo              meshInfo = MeshInfo::meshSymmetricCuboid( Point3D( {0, 0, 0} ), Point3D( {1, 1, 1} ), 1, 1, 1 );
   const SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   const auto                  storage = std::make_shared< PrimitiveStorage >( setupStorage );

   WALBERLA_CHECK_GREATER_EQUAL( setupStorage.getNumberOfCells(),
                                 uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );

   P2Function< real_t > p2Src( "x", storage, level, level );
   P2Function< real_t > p2Dst( "x", storage, level, level );

   std::function< real_t( const hhg::Point3D& ) > someFunction = []( const hhg::Point3D& x ) {
      return cos( M_PI * x[0] ) - sin( 2.0 * M_PI * x[1] ) + cos( 2.0 * M_PI * x[2] );
   };

   p2Src.interpolate( someFunction, level );
   p2Dst.interpolate( someFunction, level );

   P2ConstantLaplaceOperator p2Operator( storage, level, level );

   // get any tet from the storage
   const Cell& cell = *storage->getCell( storage->getCellIDs().front() );

   // gather all the pointers
   auto v2v_opr_data = cell.getData( p2Operator.getVertexToVertexOpr().getCellStencilID() )->getData( level );
   auto v2v_src_data = cell.getData( p2Src.getVertexDoFFunction().getCellDataID() )->getPointer( level );
   auto v2v_dst_data = cell.getData( p2Dst.getVertexDoFFunction().getCellDataID() )->getPointer( level );

   auto v2e_opr_data = cell.getData( p2Operator.getVertexToEdgeOpr().getCellStencilID() )->getData( level );
   auto v2e_src_data = cell.getData( p2Src.getVertexDoFFunction().getCellDataID() )->getPointer( level );
   auto v2e_dst_data = cell.getData( p2Dst.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

   auto e2v_opr_data = cell.getData( p2Operator.getEdgeToVertexOpr().getCellStencilID() )->getData( level );
   auto e2v_src_data = cell.getData( p2Src.getEdgeDoFFunction().getCellDataID() )->getPointer( level );
   auto e2v_dst_data = cell.getData( p2Dst.getVertexDoFFunction().getCellDataID() )->getPointer( level );

   auto e2e_opr_data = cell.getData( p2Operator.getEdgeToEdgeOpr().getCellStencilID() )->getData( level );
   auto e2e_src_data = cell.getData( p2Src.getEdgeDoFFunction().getCellDataID() )->getPointer( level );
   auto e2e_dst_data = cell.getData( p2Dst.getEdgeDoFFunction().getCellDataID() )->getPointer( level );

   typedef edgedof::EdgeDoFOrientation eo;
   std::map< eo, uint_t >              firstEdgeIdx;
   for ( auto e : edgedof::allEdgeDoFOrientations )
      firstEdgeIdx[e] = edgedof::macrocell::index( level, 0, 0, 0, e );

   /////////////////////
   // Apply benchmark //
   /////////////////////

   uint_t chunk = 0;
   while ( timer.total() < iterationMinTime )
   {
      timer.start();
      LIKWID_MARKER_START( "OperatorBenchmark" );
      switch ( kernelType )
      {
      case APPLY_V_TO_V_REPLACE:

         for ( uint_t i = 0; i < chunkIterations; i++ )
         {
            vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_replace(
                v2v_dst_data, v2v_src_data, static_cast< int32_t >( level ), v2v_opr_data );
            hhg::misc::dummy( v2v_src_data );
            hhg::misc::dummy( v2v_dst_data );
         }
         break;

      case APPLY_V_TO_V_ADD:
         for ( uint_t i = 0; i < chunkIterations; i++ )
         {
            vertexdof::macrocell::generated::apply_3D_macrocell_vertexdof_to_vertexdof_add(
                v2v_dst_data, v2v_src_data, static_cast< int32_t >( level ), v2v_opr_data );
            hhg::misc::dummy( v2v_src_data );
            hhg::misc::dummy( v2v_dst_data );
         }
         break;

      case APPLY_E_TO_V_REPLACE:
         for ( uint_t i = 0; i < chunkIterations; i++ )
         {
            EdgeDoFToVertexDoF::generated::apply_3D_macrocell_edgedof_to_vertexdof_replace( &e2v_src_data[firstEdgeIdx[eo::X]],
                                                                                            &e2v_src_data[firstEdgeIdx[eo::XY]],
                                                                                            &e2v_src_data[firstEdgeIdx[eo::XYZ]],
                                                                                            &e2v_src_data[firstEdgeIdx[eo::XZ]],
                                                                                            &e2v_src_data[firstEdgeIdx[eo::Y]],
                                                                                            &e2v_src_data[firstEdgeIdx[eo::YZ]],
                                                                                            &e2v_src_data[firstEdgeIdx[eo::Z]],
                                                                                            e2v_dst_data,
                                                                                            e2v_opr_data,
                                                                                            static_cast< int32_t >( level ) );
            hhg::misc::dummy( e2v_src_data );
            hhg::misc::dummy( e2v_dst_data );
         }
         break;

      case APPLY_E_TO_E_REPLACE:
         for ( uint_t i = 0; i < chunkIterations; i++ )
         {
            edgedof::macrocell::generated::apply_3D_macrocell_edgedof_to_edgedof_replace( &e2e_dst_data[firstEdgeIdx[eo::X]],
                                                                                          &e2e_dst_data[firstEdgeIdx[eo::XY]],
                                                                                          &e2e_dst_data[firstEdgeIdx[eo::XYZ]],
                                                                                          &e2e_dst_data[firstEdgeIdx[eo::XZ]],
                                                                                          &e2e_dst_data[firstEdgeIdx[eo::Y]],
                                                                                          &e2e_dst_data[firstEdgeIdx[eo::YZ]],
                                                                                          &e2e_dst_data[firstEdgeIdx[eo::Z]],
                                                                                          &e2e_src_data[firstEdgeIdx[eo::X]],
                                                                                          &e2e_src_data[firstEdgeIdx[eo::XY]],
                                                                                          &e2e_src_data[firstEdgeIdx[eo::XYZ]],
                                                                                          &e2e_src_data[firstEdgeIdx[eo::XZ]],
                                                                                          &e2e_src_data[firstEdgeIdx[eo::Y]],
                                                                                          &e2e_src_data[firstEdgeIdx[eo::YZ]],
                                                                                          &e2e_src_data[firstEdgeIdx[eo::Z]],
                                                                                          e2e_opr_data,
                                                                                          static_cast< int32_t >( level ) );
            hhg::misc::dummy( e2e_src_data );
            hhg::misc::dummy( e2e_dst_data );
         }
         break;

      default:
         WALBERLA_ABORT( "Kernel type not implemented." );
         break;
      }
      LIKWID_MARKER_STOP( "OperatorBenchmark" );
      timer.end();

      chunk++;
      WALBERLA_LOG_INFO_ON_ROOT( "Iteration: " << chunk * chunkIterations << ": " << timer.last() << "sec" );
   }

   LIKWID_MARKER_CLOSE;

   WALBERLA_LOG_INFO_ON_ROOT( "Total number of iterations: " << chunk * chunkIterations << ", total time: " << timer.total()
                                                             << "sec" );
}
} // namespace hhg
int main( int argc, char** argv )
{
   LIKWID_MARKER_INIT;

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   LIKWID_MARKER_THREADINIT;

   LIKWID_MARKER_REGISTER( "OperatorBenchmark" );

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./P2OperatorBenchmarks.prm";
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }

   const walberla::Config::BlockHandle mainConf         = cfg->getBlock( "Parameters" );
   const uint_t                        level            = mainConf.getParameter< uint_t >( "level" );
   const real_t                        iterationMinTime = mainConf.getParameter< real_t >( "iterationMinTime" );
   const std::string                   kernelType       = mainConf.getParameter< std::string >( "kernelType" );

   hhg::runBenchmark( level, iterationMinTime, hhg::strToKernelType.at( kernelType ) );
}