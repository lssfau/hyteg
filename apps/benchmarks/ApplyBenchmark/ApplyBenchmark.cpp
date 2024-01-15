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
#include "core/config/Config.h"
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"
#include "core/timing/TimingJSON.h"

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

using walberla::real_c;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

template < typename Discretization, typename Operator >
static walberla::WcTimingTree runbenchmark( const uint_t& level, const uint_t& facesPerProcess, const uint_t& flopsPerIter )
{
   const uint_t  numProc  = uint_c( walberla::mpi::MPIManager::instance()->numProcesses() );
   hyteg::MeshInfo meshInfo = hyteg::MeshInfo::meshFaceChain( numProc * facesPerProcess );

   hyteg::SetupPrimitiveStorage setupStorage( meshInfo, numProc );
   hyteg::loadbalancing::roundRobin( setupStorage );

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );
   std::shared_ptr< hyteg::PrimitiveStorage >  storage = std::make_shared< hyteg::PrimitiveStorage >( setupStorage, timingTree );

   std::function< real_t( const hyteg::Point3D& ) > exact = []( const hyteg::Point3D& xx ) {
      return std::sin( walberla::math::pi * xx[0] ) + std::cos( walberla::math::pi * xx[1] );
   };

   auto storageInfo = storage->getGlobalInfo();
   WALBERLA_LOG_PROGRESS_ON_ROOT( storageInfo )

   ///// Functions / operators / allocation /////

   Discretization src( "src", storage, level, level );
   Discretization dst( "dst", storage, level, level );

   Operator laplace( storage, level, level );

   [[maybe_unused]]const uint_t localDoFs = hyteg::numberOfLocalDoFs< typename Discretization::Tag >( *storage, level );
   [[maybe_unused]]const uint_t totalDoFs = hyteg::numberOfGlobalDoFs< typename Discretization::Tag >( *storage, level );

   src.interpolate( exact, level, hyteg::Inner );

   WALBERLA_LOG_PROGRESS( "localDoFs: " << localDoFs << " totalDoFs: " << totalDoFs )

   walberla::WcTimer timer;
   uint_t            iterations = 1;
   do
   {
      LIKWID_MARKER_RESET( "HyTeG-apply" );
      LIKWID_MARKER_START( "HyTeG-apply" );
      timer.reset();
      for ( uint_t i = 0; i < iterations; ++i )
      {
         laplace.apply( src, dst, level, hyteg::Inner );
      }
      timer.end();
      LIKWID_MARKER_STOP( "HyTeG-apply" );

      iterations *= 2;
   } while ( timer.last() < 1 );

   iterations /= 2;

   double hyteg_apply = timer.last();
   WALBERLA_LOG_PROGRESS_ON_ROOT( "HyTeG apply runtime: " << hyteg_apply )

   walberla::WcTimingTree tt  = timingTree->getReduced();
   auto                   tt2 = tt.getCopyWithRemainder();

   const uint_t globalInnerDoFs = hyteg::numberOfGlobalInnerDoFs< typename Discretization::Tag >( *storage, level );
   const real_t glups           = real_c( globalInnerDoFs * iterations ) / 1e9 / hyteg_apply;
   const real_t gflops          = real_c( globalInnerDoFs * iterations * flopsPerIter ) / 1e9 / hyteg_apply;

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%10.3e|%10.3e|%10.3e|%10.3e|%5u|%5u|%7u",
                                           hyteg_apply,
                                           glups,
                                           gflops,
                                           real_c( globalInnerDoFs ),
                                           level,
                                           numProc,
                                           facesPerProcess ) )

   return tt2;
}
int main( int argc, char* argv[] )
{
   LIKWID_MARKER_INIT;
   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   LIKWID_MARKER_THREADINIT;

   auto cfg = std::make_shared< walberla::config::Config >();
   if ( env.config() == nullptr )
   {
      auto defaultFile = "./ApplyBenchmark.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT( "No Parameter file given loading default parameter file: " << defaultFile )
      cfg->readParameterFile( defaultFile );
   }
   else
   {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf        = cfg->getBlock( "Parameters" );
   const uint_t                        level           = mainConf.getParameter< uint_t >( "level" );
   const uint_t                        facesPerProcess = mainConf.getParameter< uint_t >( "facesPerProcess" );
   const uint_t                        logLevel        = mainConf.getParameter< uint_t >( "logLevel" );
   const std::string                   discretization  = mainConf.getParameter< std::string >( "discretization" );

   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::LogLevel( logLevel ) );

   walberla::WcTimingTree tt2;

   WALBERLA_LOG_INFO_ON_ROOT( walberla::format( "%10s|%10s|%10s|%10s|%5s|%5s|%7s| Discr.: %s",
                                           "Time (s)",
                                           "GDoF/s",
                                           "GFLOP/s",
                                           " DoFs ",
                                           "Level",
                                           "Procs",
                                           "face/proc",
                                           discretization.c_str() ) )

   if ( discretization == "P1" )
   {
      tt2 = runbenchmark< hyteg::P1Function< real_t >, hyteg::P1ConstantLaplaceOperator >( level, facesPerProcess, 13 );
   }
   else if ( discretization == "P2" )
   {
      const uint_t flops           = 13 + 21 + 23 + 27;
      tt2 = runbenchmark< hyteg::P2Function< real_t >, hyteg::P2ConstantLaplaceOperator >( level, facesPerProcess, flops );
   }
   else
   {
      WALBERLA_ABORT( "Unknown discretization: " << discretization )
   }

   if ( mainConf.getParameter< bool >( "printTiming" ) )
   {
      WALBERLA_LOG_INFO_ON_ROOT( tt2 )
   }

   if ( mainConf.getParameter< bool >( "writeJSON" ) )
   {
      nlohmann::json ttjson = nlohmann::json( tt2 );
      std::ofstream  o( "ApplyBenchmarkOutput.json" );
      o << ttjson;
      o.close();
   }

   LIKWID_MARKER_CLOSE;
}
