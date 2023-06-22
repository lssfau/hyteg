/*
 * Copyright (c) 2017-2019 Christoph Schwarzmeier, Dominik Thoennes, Nils Kohl.
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
#include "core/Format.hpp"
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"
#include "core/timing/TimingJSON.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/edgedofspace/generatedKernels/apply_2D_macroface_edgedof_to_edgedof_replace.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/misc/dummy.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/EdgeDoFToVertexDoFApply.hpp"
#include "hyteg/mixedoperators/EdgeDoFToVertexDoFOperator/generatedKernels/apply_2D_macroface_edgedof_to_vertexdof_replace.hpp"
#include "hyteg/mixedoperators/VertexDoFToEdgeDoFOperator/generatedKernels/apply_2D_macroface_vertexdof_to_edgedof_replace.hpp"
#include "hyteg/p1functionspace/generatedKernels/apply_2D_macroface_vertexdof_to_vertexdof_replace.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "AbstractApply.hpp"

using walberla::real_t;
using namespace hyteg;

///global varialbes controling which benchmarks to execute
bool benchVertexToVertex = true;
bool benchEdgeToVertex   = true;
bool benchEdgeToEdge     = true;
bool benchVertextoEdge   = true;

static void registerLikwidRegion( const std::string& regionName )
{
   WALBERLA_UNUSED( regionName );
   /// register, start and stop to avoid warning in RESET
   LIKWID_MARKER_REGISTER( regionName.c_str() );
   LIKWID_MARKER_START( regionName.c_str() );
   LIKWID_MARKER_STOP( regionName.c_str() );
}

static void performBenchmark( hyteg::P2Function< double >&      src,
                              hyteg::P2Function< double >&      dst,
                              hyteg::P2ConstantLaplaceOperator& laplace,
                              const uint_t&                     level,
                              Face&                             face,
                              walberla::WcTimingTree&           timingTree,
                              const bool                        useGeneratedKernels,
                              uint_t                            startIterations )
{
   const std::string benchInfoString = "level" + ( level < 10 ? "0" + std::to_string( level ) : std::to_string( level ) ) +
                                       "-numProcs" + std::to_string( walberla::mpi::MPIManager::instance()->numProcesses() );

   double time = 0, mlups, mflops;

#ifdef LIKWID_PERFMON
   double events;
   int    nevents = 0, count;
#endif

   uint_t      iterations = startIterations;
   std::string vvname, vename, eename, evname;
   vvname = "Vertex-to-Vertex-Apply-" + benchInfoString;
   evname = "Edge-to-Vertex-Apply-" + benchInfoString;
   eename = "Edge-to-Edge-Apply-" + benchInfoString;
   vename = "Vertex-to-Edge-Apply-" + benchInfoString;

   uint_t innerIterationsVertex =
       levelinfo::num_microvertices_per_face_from_width( levelinfo::num_microvertices_per_edge( level ) - 3 );
   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( "%18s|%10s|%10s|%10s|%6s|%5s", "kernel", "Time (s)", "MLUPs", "MFLOPs", " Iter", " Level" ) )

   typedef edgedof::EdgeDoFOrientation eo;
   std::map< eo, uint_t >              firstIdx;
   for ( auto e : edgedof::faceLocalEdgeDoFOrientations )
      firstIdx[e] = edgedof::macroface::index( level, 0, 0, e );

   /// Vertex to Vertex

   if ( benchVertexToVertex )
   {
      registerLikwidRegion( vvname );
      do
      {
         timingTree.start( vvname );
         ///only works with likwid 4.3.3 and higher
         LIKWID_MARKER_RESET( vvname.c_str() );
         LIKWID_MARKER_START( vvname.c_str() );

         auto dstPtr     = face.getData( dst.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         auto srcPtr     = face.getData( src.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         auto stencilPtr = face.getData( laplace.getVertexToVertexOpr().getFaceStencilID() )->getPointer( level );

         for ( uint_t iter = 0; iter < iterations; ++iter )
         {
            if ( useGeneratedKernels )
            {
               hyteg::vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_replace(
                   dstPtr, srcPtr, stencilPtr, static_cast< int32_t >( level ) );
            }
            else
            {
               apply_2d_vertex_to_vertex( dstPtr, srcPtr, stencilPtr, level );
            }
         }
         WALBERLA_MPI_BARRIER()
         LIKWID_MARKER_STOP( vvname.c_str() );
         timingTree.stop( vvname );
         iterations *= 2;
      } while ( timingTree[vvname].last() < 0.5 );

      iterations /= 2;

#ifdef LIKWID_PERFMON
      LIKWID_MARKER_GET( vvname.c_str(), &nevents, &events, &time, &count );
      if ( time < 1e-16 )
      {
         time = timingTree[vvname].last();
      }
#else
      time = timingTree[vvname].last();
#endif

      mlups = real_t( innerIterationsVertex * iterations ) / time / 1e6;
      /// 13 Flops: 7 Mults and 6 Adds
      mflops = real_t( innerIterationsVertex * iterations * 13 ) / time / 1e6;

      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( "%18s|%10.3e|%10.3e|%10.3e|%6u|%5u", "vertex to vertex", time, mlups, mflops, iterations, level ) )
   }
   /// Edge to Vertex
   if ( benchEdgeToVertex )
   {
      registerLikwidRegion( evname );
      iterations = startIterations;
      do
      {
         timingTree.start( evname );
         ///only works with likwid 4.3.3 and higher
         LIKWID_MARKER_RESET( evname.c_str() );
         LIKWID_MARKER_START( evname.c_str() );

         auto dstPtr     = face.getData( dst.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         auto srcPtr     = face.getData( src.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
         auto stencilPtr = face.getData( laplace.getEdgeToVertexOpr().getFaceStencilID() )->getPointer( level );

         for ( uint_t i = 0; i < iterations; i++ )
         {
            hyteg::EdgeDoFToVertexDoF::generated::apply_2D_macroface_edgedof_to_vertexdof_replace(
                &srcPtr[firstIdx[eo::X]],
                &srcPtr[firstIdx[eo::XY]],
                &srcPtr[firstIdx[eo::Y]],
                stencilPtr,
                dstPtr,
                static_cast< int32_t >( level ) );
            hyteg::misc::dummy( srcPtr, dstPtr );
         }
         WALBERLA_MPI_BARRIER()
         LIKWID_MARKER_STOP( evname.c_str() );
         timingTree.stop( evname );
         iterations *= 2;
      } while ( timingTree[evname].last() < 0.5 );

      iterations /= 2;

#ifdef LIKWID_PERFMON
      LIKWID_MARKER_GET( evname.c_str(), &nevents, &events, &time, &count );
      if ( time < 1e-16 )
      {
         time = timingTree[evname].last();
      }
#else
      time = timingTree[evname].last();
#endif
      mlups = real_t( innerIterationsVertex * iterations ) / time / 1e6;
      /// 4 DoFs for each subgroup; 23 Flops: 12 Mults and 11 Adds
      mflops = real_t( innerIterationsVertex * iterations * 23 ) / time / 1e6;

      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( "%18s|%10.3e|%10.3e|%10.3e|%6u|%5u", "edge to vertex", time, mlups, mflops, iterations, level ) )
   }
   /// Edge to Edge
   if ( benchEdgeToEdge )
   {
      registerLikwidRegion( eename );
      iterations = startIterations;
      do
      {
         timingTree.start( eename );
         ///only works with likwid 4.3.3 and higher
         LIKWID_MARKER_RESET( eename.c_str() );
         LIKWID_MARKER_START( eename.c_str() );

         auto dstPtr     = face.getData( dst.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
         auto srcPtr     = face.getData( src.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
         auto stencilPtr = face.getData( laplace.getEdgeToEdgeOpr().getFaceStencilID() )->getPointer( level );

         for ( uint_t i = 0; i < iterations; i++ )
         {
            hyteg::edgedof::macroface::generated::apply_2D_macroface_edgedof_to_edgedof_replace(
                &dstPtr[firstIdx[eo::X]],
                &dstPtr[firstIdx[eo::XY]],
                &dstPtr[firstIdx[eo::Y]],
                &srcPtr[firstIdx[eo::X]],
                &srcPtr[firstIdx[eo::XY]],
                &srcPtr[firstIdx[eo::Y]],
                &stencilPtr[5],
                &stencilPtr[0],
                &stencilPtr[10],
                static_cast< int32_t >( level ) );
            hyteg::misc::dummy( srcPtr, dstPtr );
         }
         WALBERLA_MPI_BARRIER()
         LIKWID_MARKER_STOP( eename.c_str() );
         timingTree.stop( eename );
         iterations *= 2;
      } while ( timingTree[eename].last() < 0.5 );

      iterations /= 2;

#ifdef LIKWID_PERFMON
      LIKWID_MARKER_GET( eename.c_str(), &nevents, &events, &time, &count );
      if ( time < 1e-16 )
      {
         time = timingTree[eename].last();
      }
#else
      time = timingTree[eename].last();
#endif
      mlups = real_t( innerIterationsVertex * iterations ) / time / 1e6;
      /// 5 DoFs for each subgroup; 29 Flops: 15 Mults and 12 Adds
      mflops = real_t( innerIterationsVertex * iterations * 27 ) / time / 1e6;

      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( "%18s|%10.3e|%10.3e|%10.3e|%6u|%5u", "edge to edge", time, mlups, mflops, iterations, level ) )
   }
   /// Vertex to Edge
   if ( benchVertextoEdge )
   {
      registerLikwidRegion( vename );
      iterations = startIterations;
      do
      {
         timingTree.start( vename );
         ///only works with likwid 4.3.3 and higher
         LIKWID_MARKER_RESET( vename.c_str() );
         LIKWID_MARKER_START( vename.c_str() );

         auto dstPtr                      = face.getData( dst.getEdgeDoFFunction().getFaceDataID() )->getPointer( level );
         auto srcPtr                      = face.getData( src.getVertexDoFFunction().getFaceDataID() )->getPointer( level );
         auto stencilPtr                  = face.getData( laplace.getVertexToEdgeOpr().getFaceStencilID() )->getPointer( level );
         auto vertexToDiagonalEdgeStencil = &stencilPtr[4];
         auto vertexToHorizontalEdgeStencil = &stencilPtr[0];
         auto vertexToVerticalEdgeStencil   = &stencilPtr[8];

         for ( uint_t i = 0; i < iterations; i++ )
         {
            hyteg::VertexDoFToEdgeDoF::generated::apply_2D_macroface_vertexdof_to_edgedof_replace(
                &dstPtr[firstIdx[eo::X]],
                &dstPtr[firstIdx[eo::XY]],
                &dstPtr[firstIdx[eo::Y]],
                srcPtr,
                vertexToDiagonalEdgeStencil,
                vertexToHorizontalEdgeStencil,
                vertexToVerticalEdgeStencil,
                static_cast< int32_t >( level ) );
            hyteg::misc::dummy( srcPtr, dstPtr );
         }
         WALBERLA_MPI_BARRIER()
         LIKWID_MARKER_STOP( vename.c_str() );
         timingTree.stop( vename );
         iterations *= 2;
      } while ( timingTree[vename].last() < 0.5 );

      iterations /= 2;

#ifdef LIKWID_PERFMON
      LIKWID_MARKER_GET( vename.c_str(), &nevents, &events, &time, &count );
      if ( time < 1e-16 )
      {
         time = timingTree[vename].last();
      }
#else
      time = timingTree[vename].last();
#endif
      mlups = real_t( innerIterationsVertex * iterations ) / time / 1e6;
      /// 21 Flops: 12 Mults and 3 * 3 Adds
      mflops = real_t( innerIterationsVertex * iterations * 21 ) / time / 1e6;

      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( "%18s|%10.3e|%10.3e|%10.3e|%6u|%5u", "vertex to edge", time, mlups, mflops, iterations, level ) )
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
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./ApplyPerformanceAnalysis-2D-P2.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile )
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const uint_t              level               = mainConf.getParameter< uint_t >( "level" );
   const bool                useGeneratedKernels = mainConf.getParameter< bool >( "useGenKernel" );
   const uint_t              startIterations     = mainConf.getParameter< uint_t >( "startIterations" );
   benchVertexToVertex                           = mainConf.getParameter< bool >( "VertexToVertex" );
   benchEdgeToVertex                             = mainConf.getParameter< bool >( "EdgeToVertex" );
   benchEdgeToEdge                               = mainConf.getParameter< bool >( "EdgeToEdge" );
   benchVertextoEdge                             = mainConf.getParameter< bool >( "VertextoEdge" );

   hyteg::MeshInfo meshInfo = hyteg::MeshInfo::meshFaceChain( uint_c( walberla::MPIManager::instance()->numProcesses() ) );
   hyteg::SetupPrimitiveStorage setupStorage( meshInfo,
                                              walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   hyteg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< walberla::WcTimingTree >  timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hyteg::PrimitiveStorage > storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) {
      return std::sin( walberla::math::pi * xx[0] ) + std::cos( walberla::math::pi * xx[1] );
   };

   ///// Functions / operators / allocation /////

   hyteg::P2Function< double > src( "src", storage, level, level );
   hyteg::P2Function< double > dst( "dst", storage, level, level );

   hyteg::P2ConstantLaplaceOperator laplace( storage, level, level );

   src.interpolate( exact, level, hyteg::Inner );

   ///// Apply benchmarks /////

   std::vector< PrimitiveID > macroFaces = storage->getFaceIDs();
   WALBERLA_CHECK_EQUAL( macroFaces.size(), 1 )
   auto face = storage->getFace( macroFaces.front() );

   performBenchmark( src, dst, laplace, level, *face, wcTimingTreeApp, useGeneratedKernels, startIterations );

   if ( mainConf.getParameter< bool >( "printTiming" ) )
   {
      auto wcTPReduced = wcTimingTreeApp.getReduced();
      WALBERLA_LOG_INFO_ON_ROOT( wcTPReduced )

      walberla::WcTimingTree tt  = timingTree->getReduced();
      auto                   tt2 = tt.getCopyWithRemainder();
      WALBERLA_LOG_INFO_ON_ROOT( tt2 )

      nlohmann::json ttJson;
      walberla::timing::to_json( ttJson, wcTPReduced );
      std::ofstream jsonOutput;
      jsonOutput.open( "ApplyPerformanceAnalysis-2D-P2.json" );
      jsonOutput << ttJson.dump( 4 );
      jsonOutput.close();
   }

   if ( mainConf.getParameter< bool >( "VTKOutput" ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Writing VTK output" )
      hyteg::VTKOutput vtkOutput( "./output", "ApplyPerformanceAnalysis-2D-P2.cpp", storage );
      vtkOutput.add( src );
      vtkOutput.add( dst );
      vtkOutput.write( level );
   }

   LIKWID_MARKER_CLOSE;
}
