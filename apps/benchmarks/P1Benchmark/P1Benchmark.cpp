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

#include "hyteg/LikwidWrapper.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroCell.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroEdge.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroFace.hpp"
#include "hyteg/p1functionspace/VertexDoFMacroVertex.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"
#include "constant_stencil_operator/P1generatedKernels/apply_2D_macroface_vertexdof_to_vertexdof_replace.hpp"

using walberla::real_c;
using walberla::real_t;
using namespace hyteg;

void kernel( double* fd_p1FaceDst, double* fd_p1FaceSrc, double* fd_p1FaceStencil, uint_t level )
{
   for( int ctr_2 = 1; ctr_2 < ( 1 << level ) - 1; ctr_2 += 1 )
      for( int ctr_1 = 1; ctr_1 < -ctr_2 + ( 1 << level ); ctr_1 += 1 )
      {
         fd_p1FaceDst[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 )] =
             fd_p1FaceSrc[ctr_1 + ( ctr_2 + 1 ) * ( ( 1 << level ) + 2 ) - ( ( ctr_2 + 1 ) * ( ctr_2 + 2 ) / 2 ) - 1] *
                 fd_p1FaceStencil[5] +
             fd_p1FaceSrc[ctr_1 + ( ctr_2 + 1 ) * ( ( 1 << level ) + 2 ) - ( ( ctr_2 + 1 ) * ( ctr_2 + 2 ) / 2 )] *
                 fd_p1FaceStencil[6] +
             fd_p1FaceSrc[ctr_1 + ( ctr_2 - 1 ) * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 - 1 ) / 2 ) + 1] *
                 fd_p1FaceStencil[1] +
             fd_p1FaceSrc[ctr_1 + ( ctr_2 - 1 ) * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 - 1 ) / 2 )] * fd_p1FaceStencil[0] +
             fd_p1FaceSrc[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 ) + 1] * fd_p1FaceStencil[4] +
             fd_p1FaceSrc[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 ) - 1] * fd_p1FaceStencil[2] +
             fd_p1FaceSrc[ctr_1 + ctr_2 * ( ( 1 << level ) + 2 ) - ( ctr_2 * ( ctr_2 + 1 ) / 2 )] * fd_p1FaceStencil[3];
      }
}
/*
 * This benchmark meassures the time for several vertexDof/P1 functions on a macro face
 */
int main( int argc, char** argv )
{
   LIKWID_MARKER_INIT;

   walberla::Environment env( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();
   walberla::logging::Logging::instance()->setStreamLogLevel( walberla::logging::Logging::DETAIL );

   //check if a config was given on command line or load default file otherwise
   auto cfg = std::make_shared<walberla::config::Config>();
   if( env.config() == nullptr ) {
      auto defaultFile = "./P1Benchmark.prm";
      WALBERLA_LOG_PROGRESS_ON_ROOT("No Parameter file given loading default parameter file: " << defaultFile);
      cfg->readParameterFile( defaultFile );
   } else {
      cfg = env.config();
   }
   const walberla::Config::BlockHandle mainConf = cfg->getBlock( "Parameters" );

   const uint_t level = mainConf.getParameter< uint_t >( "level" );

   WALBERLA_LOG_INFO("level: " << level);
   const uint_t totalInnerPoints =
       hyteg::indexing::layout::linearMacroFaceSize( levelinfo::num_microvertices_per_edge( level ) - 3 );
   WALBERLA_LOG_INFO_ON_ROOT("Total inner points: " << totalInnerPoints);

   LIKWID_MARKER_THREADINIT;

   MeshInfo                            meshInfo = MeshInfo::meshUnitSquare( 0 );
   SetupPrimitiveStorage               setupStorage( meshInfo, 1 );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   auto src = std::make_shared< hyteg::P1Function< real_t > >( "src", storage, level, level );
   auto dst = std::make_shared< hyteg::P1Function< real_t > >( "dst", storage, level, level );

   hyteg::P1ConstantLaplaceOperator M( storage, level, level );

   std::shared_ptr< Face > face = storage->getFaces().begin().operator*().second;

   real_t* stencilPtr = face->getData( M.getFaceStencilID() )->getPointer( level );

   for( uint_t i = 0; i < face->getData( M.getFaceStencilID() )->getSize( level ); ++i )
   {
      stencilPtr[i] = 1. / real_c( i % 10 + 2 );
   }

   for( uint_t i = 0; i < face->getData( M.getFaceStencilID() )->getSize( level ); ++i )
   {
      WALBERLA_LOG_DETAIL_ON_ROOT( i << ": " << stencilPtr[i] );
   }

   std::function< real_t( const hyteg::Point3D& ) > exactFunc = [&]( const hyteg::Point3D& point ) {
      return sqrt( point[0] * point[0] + point[1] * point[1] );
   };
   std::function< real_t( const Point3D&, const std::vector< real_t >& ) > exactExtended =
         [&exactFunc]( const hyteg::Point3D& x, const std::vector< real_t >& ) { return exactFunc( x ); };

   //P1Function< real_t > x("x", storage, level, level);
   src->interpolate( exactFunc, level );

   walberla::WcTimer timer;

   LIKWID_MARKER_START( "interpolate" );
   timer.reset();
   vertexdof::macroface::interpolate< real_t >(
         level, *face, src->getFaceDataID(), {}, exactExtended );
   timer.end();
   LIKWID_MARKER_STOP( "interpolate" );
   WALBERLA_LOG_INFO_ON_ROOT( "interpolate runtime: " << timer.last() );

   LIKWID_MARKER_START( "apply" );
   if( hyteg::globalDefines::useGeneratedKernels ){
      WALBERLA_LOG_INFO_ON_ROOT("Using generated 2D apply kernel");
      real_t* opr_data = face->getData( M.getFaceStencilID() )->getPointer( level );
      real_t* src_data      = face->getData( src->getFaceDataID() )->getPointer( level );
      real_t* dst_data      = face->getData( dst->getFaceDataID() )->getPointer( level );
      timer.reset();
      vertexdof::macroface::generated::apply_2D_macroface_vertexdof_to_vertexdof_replace( dst_data, src_data, opr_data, static_cast< int32_t >( level ) );
   } else {
      vertexdof::macroface::apply< real_t >(
            level, *face, M.getFaceStencilID(), src->getFaceDataID(), dst->getFaceDataID(), Replace );
   }
   timer.end();
   LIKWID_MARKER_STOP( "apply" );
   WALBERLA_LOG_INFO_ON_ROOT( "apply runtime: " << timer.last() );

   real_t check1 = vertexdof::macroface::dot< real_t >( level, *face, dst->getFaceDataID(), dst->getFaceDataID() );

   LIKWID_MARKER_START( "assign" );
   timer.reset();
   vertexdof::macroface::assign< real_t >( level, *face, {13}, {src->getFaceDataID()}, dst->getFaceDataID() );
   timer.end();
   LIKWID_MARKER_STOP( "assign" );
   WALBERLA_LOG_INFO_ON_ROOT( "assign runtime: " << timer.last() );

   real_t* opr_data = face->getData( M.getFaceStencilID() )->getPointer( level );
   real_t* srcPtr   = face->getData( src->getFaceDataID() )->getPointer( level );
   real_t* dstPtr   = face->getData( dst->getFaceDataID() )->getPointer( level );

   LIKWID_MARKER_START( "apply_gen" );
   timer.reset();
   kernel( dstPtr, srcPtr, opr_data, level );
   timer.end();
   LIKWID_MARKER_STOP( "apply_gen" );
   // WALBERLA_LOG_INFO_ON_ROOT( "apply_gen runtime: " << timer.last() );

   /// do something with the result to prevent the compiler from removing all the computations
   real_t check2 = vertexdof::macroface::dot< real_t >( level, *face, dst->getFaceDataID(), dst->getFaceDataID() );
   WALBERLA_LOG_INFO_ON_ROOT( check2 );
   WALBERLA_CHECK_FLOAT_EQUAL( check1, check2 );

   LIKWID_MARKER_CLOSE;
}
