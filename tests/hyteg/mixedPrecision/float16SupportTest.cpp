/*
* Copyright (c) 2023 Michael Zikeli.
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

#include "core/DataTypes.h"
#include "core/Environment.h"
#include "core/logging/Logging.h"

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

namespace hyteg::simple_Float16_test {
using namespace hyteg;
using walberla::floatIsEqual;
using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

// === Choosing Accuracy ===
//+++ Precision : fp16 +++
using walberla::float16;
using walberla::float32;
using walberla::float64;
using dst_t                         = float16;
using src_t                         = real_t;
constexpr real_t     precisionLimit = walberla::float16( 1e-3 );
const std::string    precisionType  = "float16";
constexpr const auto maxLevel       = uint_t( 3 );

void vertexDoFFct_test( const std::shared_ptr< PrimitiveStorage >& storage )
{
   vertexdof::VertexDoFFunction< src_t > fpSrc( "fpSrc", storage, maxLevel, maxLevel );
   vertexdof::VertexDoFFunction< src_t > oneFunction( "oneFunction", storage, maxLevel, maxLevel );
   vertexdof::VertexDoFFunction< src_t > err( "err", storage, maxLevel, maxLevel );
   vertexdof::VertexDoFFunction< src_t > fpDst_extend( "fpDst_extend", storage, maxLevel, maxLevel );
   vertexdof::VertexDoFFunction< dst_t > fpDst( "fpDst", storage, maxLevel, maxLevel );

   auto oneSrc = []( const Point3D& ) { return (src_t) ( 1.0 ); };
   auto oneDst = []( const Point3D& ) { return (dst_t) ( 1.0 ); };

   oneFunction.interpolate( oneSrc, maxLevel, DoFType::All );
   fpSrc.interpolate( oneSrc, maxLevel, DoFType::All );
   fpDst.interpolate( oneDst, maxLevel, DoFType::All );
   fpDst_extend.copyFrom( fpSrc, maxLevel );

   err.assign( { 1, -1 }, { fpSrc, fpDst_extend }, maxLevel );
   const auto numPointsErr = oneFunction.dotGlobal( oneFunction, maxLevel, DoFType::All );

   const auto l2error = std::sqrt( err.dotGlobal( err, maxLevel, DoFType::All ) / numPointsErr );

   WALBERLA_CHECK_LESS( l2error, precisionLimit );

} // simple_Float16_test:vertexDoFFct_test(storage)

void p1Fct_test( const std::shared_ptr< PrimitiveStorage >& storage )
{
   P1Function< src_t > fpSrc( "fpSrc", storage, maxLevel, maxLevel );
   P1Function< src_t > oneFunction( "oneFunction", storage, maxLevel, maxLevel );
   P1Function< src_t > err( "err", storage, maxLevel, maxLevel );
   P1Function< src_t > fpDst_extend( "fpDst_extend", storage, maxLevel, maxLevel );
   P1Function< dst_t > fpDst( "fpDst", storage, maxLevel, maxLevel );

   auto oneSrc = []( const Point3D& ) { return (src_t) ( 1.0 ); };
   auto oneDst = []( const Point3D& ) { return (dst_t) ( 1.0 ); };

   oneFunction.interpolate( oneSrc, maxLevel, DoFType::All );
   fpSrc.interpolate( oneSrc, maxLevel, DoFType::All );
   fpDst.interpolate( oneDst, maxLevel, DoFType::All );
   fpDst_extend.copyFrom( fpSrc, maxLevel );

   err.assign( { 1, -1 }, { fpSrc, fpDst_extend }, maxLevel );
   const auto numPointsErr = oneFunction.dotGlobal( oneFunction, maxLevel, DoFType::All );

   const auto l2error = std::sqrt( err.dotGlobal( err, maxLevel, DoFType::All ) / numPointsErr );

   WALBERLA_CHECK_LESS( l2error, precisionLimit );

} // simple_Float16_test::p1Fct_test(storage)

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::INFO );
   walberla::MPIManager::instance()->useWorldComm();

   WALBERLA_LOG_INFO_ON_ROOT( " This run is executed with " << precisionType );
   WALBERLA_LOG_INFO_ON_ROOT( " machine precision limit is " << precisionLimit );
   const std::string stringLine( 125, '=' );
   WALBERLA_LOG_INFO_ON_ROOT( stringLine );

   // +++ Set Parameters +++
   const std::string meshFile = prependHyTeGMeshDir( "2D/quad_2el.msh" );

   // +++Primitive Storage+++
   const auto            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   const auto storage = std::make_shared< PrimitiveStorage >( setupStorage );

   WALBERLA_LOG_INFO_ON_ROOT( " Start a test with a VertexDoFFunctions." );
   vertexDoFFct_test( storage );

   WALBERLA_LOG_INFO_ON_ROOT( " Start a test with a P1Function." );
   p1Fct_test( storage );

   return 0;
} // simple_Float16_test::main()

} // namespace hyteg::simple_Float16_test

int main( int argc, char** argv )
{
   hyteg::simple_Float16_test::main( argc, argv );

   return EXIT_SUCCESS;
} // main()
