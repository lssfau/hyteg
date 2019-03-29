
#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/mpi/MPIManager.h"
#include "core/timing/TimingJSON.h"
#include "core/math/Constants.h"

#include "tinyhhg_core/LikwidWrapper.hpp"
#include "tinyhhg_core/VTKWriter.hpp"
#include "tinyhhg_core/edgedofspace/generatedKernels/GeneratedKernelsEdgeToEdgeMacroFace2D.hpp"
#include "tinyhhg_core/mesh/MeshInfo.hpp"
#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFApply.hpp"
#include "tinyhhg_core/mixedoperators/EdgeDoFToVertexDoFOperator/generatedKernels/GeneratedKernelsEdgeToVertexMacroFace2D.hpp"
#include "tinyhhg_core/mixedoperators/VertexDoFToEdgeDoFOperator/generatedKernels/GeneratedKernelsVertexToEdgeMacroFace2D.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp"
#include "tinyhhg_core/p1functionspace/generatedKernels/GeneratedKernelsVertexToVertexMacroFace2D.hpp"
#include "tinyhhg_core/p2functionspace/P2ConstantOperator.hpp"
#include "tinyhhg_core/p2functionspace/P2Function.hpp"
#include "tinyhhg_core/primitivestorage/SetupPrimitiveStorage.hpp"
#include "tinyhhg_core/primitivestorage/loadbalancing/SimpleBalancer.hpp"

const int USE_GENERATED_KERNELS = 1;

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
                                       "sampleSize" + std::to_string( sampleSize ) +
                                       "numProcs" + std::to_string(walberla::mpi::MPIManager::instance()->numProcesses());

   std::string name;
   /// Vertex to Vertex
   for( uint_t i = 0; i < sampleSize; i++ )
   {
      name = "Vertex-to-Vertex-Apply-" + benchInfoString;
      timingTree.start( name );
      LIKWID_MARKER_START( name.c_str() );
      if( USE_GENERATED_KERNELS )
      {
         auto dstPtr     = face.getData( dst.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         auto srcPtr     = face.getData( src.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         auto stencilPtr = face.getData( laplace.getVertexToVertexOpr().getFaceStencilID() )->getPointer( level );
         hhg::vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_replace( dstPtr, srcPtr, stencilPtr, static_cast< int64_t  >( level ) );
      } else
      {
         hhg::vertexdof::macroface::apply( level,
                                           face,
                                           laplace.getVertexToVertexOpr().getFaceStencilID(),
                                           src.getVertexDoFFunction().getFaceDataID(),
                                           dst.getVertexDoFFunction().getFaceDataID(),
                                           hhg::Replace );
      }
      LIKWID_MARKER_STOP( name.c_str() );
      timingTree.stop( name );
   }

   /// Edge to Vertex
   for( uint_t i = 0; i < sampleSize; i++ )
   {
      name = "Edge-to-Vertex-Apply-" + benchInfoString;
      timingTree.start( name );
      LIKWID_MARKER_START( name.c_str() );
      if( USE_GENERATED_KERNELS )
      {
         auto dstPtr     = face.getData( dst.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         auto srcPtr     = face.getData( src.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
         auto stencilPtr = face.getData( laplace.getEdgeToVertexOpr().getFaceStencilID() )->getPointer( level );
         hhg::EdgeDoFToVertexDoF::generated::apply_2D_macroface_edgedof_to_vertexdof_replace( srcPtr, stencilPtr, dstPtr, static_cast< int64_t  >( level ) );
      } else
      {
         hhg::EdgeDoFToVertexDoF::applyFace( level,
                                             face,
                                             laplace.getEdgeToVertexOpr().getFaceStencilID(),
                                             src.getEdgeDoFFunction().getFaceDataID(),
                                             dst.getVertexDoFFunction().getFaceDataID(),
                                             hhg::Replace );
      }
      LIKWID_MARKER_STOP( name.c_str() );
      timingTree.stop( name );
   }

   /// Edge to Edge
   for( uint_t i = 0; i < sampleSize; i++ )
   {
      name = "Edge-to-Edge-Apply-" + benchInfoString;
      timingTree.start( name );
      LIKWID_MARKER_START( name.c_str() );
      if( USE_GENERATED_KERNELS )
      {
         auto dstPtr     = face.getData( dst.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
         auto srcPtr     = face.getData( src.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
         auto stencilPtr = face.getData( laplace.getEdgeToEdgeOpr().getFaceStencilID() )->getPointer( level );
         hhg::edgedof::macroface::generated::apply_2D_macroface_edgedof_to_edgedof_replace( dstPtr, srcPtr, &stencilPtr[5], &stencilPtr[0], &stencilPtr[10], static_cast< int64_t  >( level ) );
      } else
      {
         hhg::edgedof::macroface::apply( level,
                                         face,
                                         laplace.getEdgeToEdgeOpr().getFaceStencilID(),
                                         src.getEdgeDoFFunction().getFaceDataID(),
                                         dst.getEdgeDoFFunction().getFaceDataID(),
                                         hhg::Replace );
      }
      LIKWID_MARKER_STOP( name.c_str() );
      timingTree.stop( name );
   }

   /// Vertex to Edge
   for( uint_t i = 0; i < sampleSize; i++ )
   {
      name = "Vertex-to-Edge-Apply-" + benchInfoString;
      timingTree.start( name );
      LIKWID_MARKER_START( name.c_str() );
      if( USE_GENERATED_KERNELS )
      {
         auto dstPtr     = face.getData( dst.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
         auto srcPtr     = face.getData( src.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         auto stencilPtr = face.getData( laplace.getVertexToEdgeOpr().getFaceStencilID() )->getPointer( level );
         auto vertexToDiagonalEdgeStencil   = &stencilPtr[4];
         auto vertexToHorizontalEdgeStencil = &stencilPtr[0];
         auto vertexToVerticalEdgeStencil   = &stencilPtr[8];
         hhg::VertexDoFToEdgeDoF::generated::apply_2D_macroface_vertexdof_to_edgedof_replace( dstPtr, srcPtr, vertexToDiagonalEdgeStencil, vertexToHorizontalEdgeStencil, vertexToVerticalEdgeStencil, static_cast< int64_t  >( level ) );
      } else
      {
         hhg::VertexDoFToEdgeDoF::applyFace( level,
                                             face,
                                             laplace.getEdgeToVertexOpr().getFaceStencilID(),
                                             src.getVertexDoFFunction().getFaceDataID(),
                                             dst.getEdgeDoFFunction().getFaceDataID(),
                                             hhg::Replace );
      }
      LIKWID_MARKER_STOP( name.c_str() );
      timingTree.stop( name );
   }
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

   hhg::MeshInfo meshInfo = hhg::MeshInfo::meshFaceChain(  uint_c( walberla::MPIManager::instance()->numProcesses() ) );
   hhg::SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hhg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hhg::PrimitiveStorage > storage = std::make_shared< hhg::PrimitiveStorage >( setupStorage, timingTree );

   std::function< real_t( const hhg::Point3D& ) > exact = []( const hhg::Point3D& xx ) {
      return std::sin( walberla::math::M_PI * xx[0] ) + std::cos( walberla::math::M_PI * xx[1] );
   };

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
      hhg::VTKOutput vtkOutput("./output", "ApplyPerformanceAnalysis-2D-P2.cpp", storage);
      vtkOutput.add( src );
      vtkOutput.add( dst );
      vtkOutput.write( level );
   }

   LIKWID_MARKER_CLOSE;
}
