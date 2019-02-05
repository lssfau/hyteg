#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"
#include "core/timing/TimingJSON.h"

#include "tinyhhg_core/Format.hpp"
#include "tinyhhg_core/LikwidWrapper.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/edgedofspace/generatedKernels/GeneratedKernelsEdgeToEdgeMacroFace2D.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/misc/dummy.hpp"
#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFApply.hpp"
#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/generatedKernels/GeneratedKernelsEdgeToVertexMacroFace2D.hpp"
#include "tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/generatedKernels/GeneratedKernelsVertexToEdgeMacroFace2D.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp"
#include "tinyhhg_core/p1functionspace/generatedKernels/GeneratedKernelsVertexToVertexMacroFace2D.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_t;
using namespace hhg;

static void performBenchmark( hhg::P2Function< double >&      src,
                              hhg::P2Function< double >&      dst,
                              hhg::P2ConstantLaplaceOperator& laplace,
                              const uint_t&                   level,
                              Face&                           face,
                              const uint_t&                   sampleSize,
                              walberla::WcTimingTree&         timingTree )
{
   const std::string benchInfoString = "level" + ( level < 10 ? "0" + std::to_string( level ) : std::to_string( level ) ) + "-" +
                                       "sampleSize" + std::to_string( sampleSize ) + "numProcs" +
                                       std::to_string( walberla::mpi::MPIManager::instance()->numProcesses() );

   double      time, mlups, mflops;
   std::string vvname, vename, eename, evname;
   vvname = "Vertex-to-Vertex-Apply-" + benchInfoString;
   evname = "Edge-to-Vertex-Apply-" + benchInfoString;
   eename = "Edge-to-Edge-Apply-" + benchInfoString;
   vename = "Vertex-to-Edge-Apply-" + benchInfoString;

   LIKWID_MARKER_REGISTER( vvname.c_str() );
   LIKWID_MARKER_REGISTER( evname.c_str() );
   LIKWID_MARKER_REGISTER( eename.c_str() );
   LIKWID_MARKER_REGISTER( vename.c_str() );

   uint_t innerIterationsVertex =
       levelinfo::num_microvertices_per_face_from_width( levelinfo::num_microvertices_per_edge( level ) - 2 );
   WALBERLA_LOG_INFO_ON_ROOT( hhg::format( "%18s|%10s|%10s|%10s", "kernel", "Time (s)", "MLUPs", "MFLOPs" ) );

   /// Vertex to Vertex
   timingTree.start( vvname );
   LIKWID_MARKER_START( vvname.c_str() );
   for( uint_t i = 0; i < sampleSize; i++ )
   {
      auto dstPtr     = face.getData( dst.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
      auto srcPtr     = face.getData( src.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
      auto stencilPtr = face.getData( laplace.getVertexToVertexOpr().getFaceStencilID() )->getPointer( level );
      hhg::vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_replace(
          dstPtr, srcPtr, stencilPtr, static_cast< int64_t >( level ) );
      hhg::misc::dummy( srcPtr, dstPtr );
   }
   LIKWID_MARKER_STOP( vvname.c_str() );
   timingTree.stop( vvname );

   time  = timingTree[vvname].average();
   mlups = real_t( innerIterationsVertex * sampleSize ) / time / 1e6;
   /// 13 Flops: 7 Mults and 6 Adds
   mflops = real_t( innerIterationsVertex * sampleSize * 13 ) / time / 1e6;

   WALBERLA_LOG_INFO_ON_ROOT( hhg::format( "%18s|%10.3e|%10.3e|%10.3e", "vertex to vertex", time, mlups, mflops ) );

   /// Edge to Vertex

   timingTree.start( evname );
   LIKWID_MARKER_START( evname.c_str() );

   for( uint_t i = 0; i < sampleSize; i++ )
   {
      auto dstPtr     = face.getData( dst.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
      auto srcPtr     = face.getData( src.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
      auto stencilPtr = face.getData( laplace.getEdgeToVertexOpr().getFaceStencilID() )->getPointer( level );
      hhg::EdgeDoFToVertexDoF::generated::apply_2D_macroface_edgedof_to_vertexdof_replace(
          srcPtr, stencilPtr, dstPtr, static_cast< int64_t >( level ) );
      hhg::misc::dummy( srcPtr, dstPtr );
   }
   LIKWID_MARKER_STOP( evname.c_str() );
   timingTree.stop( evname );

   /// Edge to Edge

   timingTree.start( eename );
   LIKWID_MARKER_START( eename.c_str() );

   for( uint_t i = 0; i < sampleSize; i++ )
   {
      auto dstPtr     = face.getData( dst.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
      auto srcPtr     = face.getData( src.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
      auto stencilPtr = face.getData( laplace.getEdgeToEdgeOpr().getFaceStencilID() )->getPointer( level );
      hhg::edgedof::macroface::generated::apply_2D_macroface_edgedof_to_edgedof_replace(
          dstPtr, srcPtr, stencilPtr, static_cast< int64_t >( level ) );
      hhg::misc::dummy( srcPtr, dstPtr );
   }

   LIKWID_MARKER_STOP( eename.c_str() );
   timingTree.stop( eename );

   /// Vertex to Edge

   LIKWID_MARKER_START( vename.c_str() );
   timingTree.start( vename );

   for( uint_t i = 0; i < sampleSize; i++ )
   {
      auto dstPtr                        = face.getData( dst.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
      auto srcPtr                        = face.getData( src.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
      auto stencilPtr                    = face.getData( laplace.getVertexToEdgeOpr().getFaceStencilID() )->getPointer( level );
      auto vertexToDiagonalEdgeStencil   = &stencilPtr[4];
      auto vertexToHorizontalEdgeStencil = &stencilPtr[0];
      auto vertexToVerticalEdgeStencil   = &stencilPtr[8];
      hhg::VertexDoFToEdgeDoF::generated::apply_2D_macroface_vertexdof_to_edgedof_replace( dstPtr,
                                                                                           srcPtr,
                                                                                           vertexToDiagonalEdgeStencil,
                                                                                           vertexToHorizontalEdgeStencil,
                                                                                           vertexToVerticalEdgeStencil,
                                                                                           static_cast< int64_t >( level ) );
      hhg::misc::dummy( srcPtr, dstPtr );
   }
   LIKWID_MARKER_STOP( vename.c_str() );
   timingTree.stop( vename );
}

int main( int argc, char* argv[] )
{
   LIKWID_MARKER_INIT;

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   walberla::WcTimingTree wcTimingTreeApp;

   LIKWID_MARKER_THREADINIT;

   ///// Parameters /////

   auto cfg = std::make_shared< walberla::config::Config >();
   if( env.config() == nullptr )
   {
      auto defaultFile = "./ApplyPerformanceAnalysis-2D-P2.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile );
      cfg->readParameterFile( defaultFile );
   } else
   {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const uint_t level = mainConf.getParameter< uint_t >( "level" );

   hhg::MeshInfo meshInfo = hhg::MeshInfo::meshFaceChain( uint_c( walberla::MPIManager::instance()->numProcesses() ) );
   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hhg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage );

   std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& xx ) {
      return std::sin( walberla::math::PI * xx[0] ) + std::cos( walberla::math::PI * xx[1] );
   };

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   storage->setTimingTree( timingTree );

   ///// Functions / operators / allocation /////

   hhg::P2Function< double > src( "src", storage, level, level );
   hhg::P2Function< double > dst( "dst", storage, level, level );

   hhg::P2ConstantLaplaceOperator laplace( storage, level, level );

   src.interpolate( exact, level, hhg::Inner );

   ///// Apply benchmarks /////

   std::vector< PrimitiveID > macroFaces = storage->getFaceIDs();
   WALBERLA_CHECK_EQUAL( macroFaces.size(), 1 );
   auto face = storage->getFace( macroFaces.front() );

   WALBERLA_CHECK_LESS_EQUAL( level, 14 );
   uint_t sampleSize = ( 1u << ( 14u - level ) );
   performBenchmark( src, dst, laplace, level, *face, sampleSize, wcTimingTreeApp );

   if( mainConf.getParameter< bool >( "printTiming" ) )
   {
      auto wcTPReduced = wcTimingTreeApp.getReduced();
      WALBERLA_LOG_INFO_ON_ROOT( wcTPReduced );

      walberla::WcTimingTree tt  = timingTree->getReduced();
      auto                   tt2 = tt.getCopyWithRemainder();
      WALBERLA_LOG_INFO_ON_ROOT( tt2 );

      nlohmann::json ttJson;
      walberla::timing::to_json( ttJson, wcTPReduced );
      std::ofstream jsonOutput;
      jsonOutput.open( "ApplyPerformanceAnalysis-2D-P2.json" );
      jsonOutput << ttJson.dump( 4 );
      jsonOutput.close();
   }

   if( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Writing VTK output" );
      hhg::VTKOutput vtkOutput( "./output", "ApplyPerformanceAnalysis-2D-P2.cpp", storage );
      vtkOutput.add( src );
      vtkOutput.add( dst );
      vtkOutput.write( level );
   }

   LIKWID_MARKER_CLOSE;
}
