/*
 * Copyright (c) 2017-2020 Dominik Thoennes, Nils Kohl.
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
#include "core/Format.hpp"
#include "core/config/Config.h"
#include "core/math/Random.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"

#include "constant_stencil_operator/P1ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;

using namespace hyteg;

/// comparing constant operator and elementwise operator Jacobi
static void test( const std::string& meshFile, const uint_t& level, const uint_t& maxiter )
{
   const bool writeVTK    = false;
   const bool printTiming = true;

   std::shared_ptr< walberla::WcTimingTree > timingTree( new walberla::WcTimingTree() );

   const auto            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage, timingTree );

   P1ConstantLaplaceOperator    Lconst( storage, level, level );
   P1ElementwiseLaplaceOperator Lelem( storage, level, level );

   Lconst.computeInverseDiagonalOperatorValues();
   Lelem.computeInverseDiagonalOperatorValues();

   P1Function< real_t > residuum( "residuum", storage, level, level );
   P1Function< real_t > rhs( "rhs", storage, level, level );
   P1Function< real_t > p2functionConst( "p1FunctionConst", storage, level, level );
   P1Function< real_t > p2functionElem( "p1FunctionElem", storage, level, level );
   P1Function< real_t > Lu( "Lu", storage, level, level );
   P1Function< real_t > p2Exact( "p1Exact", storage, level, level );
   P1Function< real_t > error( "error", storage, level, level );
   P1Function< real_t > helperFun( "helperFun", storage, level, level );

   VTKOutput vtkOutput( "../../output", "Jacobi_P1", storage );
   vtkOutput.add( p2functionConst );
   vtkOutput.add( p2functionElem );
   vtkOutput.add( p2Exact );
   vtkOutput.add( rhs );
   vtkOutput.add( residuum );
   vtkOutput.add( error );
   vtkOutput.add( helperFun );

   std::function< real_t( const Point3D& ) > exactFunction = []( const Point3D& x ) { return sin( x[0] ) * sinh( x[1] ); };
   std::function< real_t( const Point3D& ) > ones          = []( const Point3D& ) { return 1.0; };
   walberla::math::seedRandomGenerator( 0 );

   p2functionConst.interpolate( exactFunction, level, DirichletBoundary );
   p2functionElem.interpolate( exactFunction, level, DirichletBoundary );
   p2Exact.interpolate( exactFunction, level );

   real_t begin_res, abs_res_old, rel_res = 0, abs_res = 0;

   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( "%6s|%10s|%10s|%10s|%10s", "iter", "Jacobi diff", "abs_res", "rel_res", "conv" ) );

   Lconst.apply( p2functionConst, Lu, level, Inner );
   residuum.assign( { 1.0, -1.0 }, { rhs, Lu }, level, Inner );
   begin_res   = std::sqrt( residuum.dotGlobal( residuum, level, Inner ) );
   abs_res_old = begin_res;

   WALBERLA_LOG_INFO_ON_ROOT(
       walberla::format( "%6d|%10.3e|%10.3e|%10.3e|%10.3e", 0, 0, begin_res, rel_res, begin_res / abs_res_old ) )
   walberla::WcTimer timer;

   for ( uint_t i = 0; i < maxiter; ++i )
   {
      if ( writeVTK )
      {
         vtkOutput.write( level, i );
      }

      helperFun.assign( { 1.0 }, { p2functionConst }, level, All );
      Lconst.smooth_jac( p2functionConst, rhs, helperFun, 1.0, level, Inner );
      helperFun.assign( { 1.0 }, { p2functionElem }, level, All );
      Lelem.smooth_jac( p2functionElem, rhs, helperFun, 1.0, level, Inner );

      error.assign( { 1.0, -1.0 }, { p2functionConst, p2functionElem }, level, All );
      auto jacobiDiff = std::sqrt( error.dotGlobal( error, level, Inner ) );

      bool dp = std::is_same< real_t, double >();
      WALBERLA_CHECK_LESS( jacobiDiff, dp ? 1e-13 : 8e-6, "ElemOp Jacobi produces different result than ConstOp Jacobi" );

      Lconst.apply( p2functionConst, Lu, level, Inner );
      residuum.assign( { 1.0, -1.0 }, { rhs, Lu }, level, Inner );
      abs_res = std::sqrt( residuum.dotGlobal( residuum, level, Inner ) );
      rel_res = abs_res / begin_res;
      WALBERLA_LOG_INFO_ON_ROOT(
          walberla::format( "%6d|%10.3e|%10.3e|%10.3e|%10.3e", i + 1, jacobiDiff, abs_res, rel_res, abs_res / abs_res_old ) )
      WALBERLA_CHECK_LESS( abs_res, abs_res_old );
      abs_res_old = abs_res;
   }
   timer.end();

   WALBERLA_LOG_INFO_ON_ROOT( "time was: " << timer.last() );
   error.assign( { 1.0, -1.0 }, { p2functionConst, p2Exact }, level );

   helperFun.interpolate( ones, level );
   real_t npoints = helperFun.dotGlobal( helperFun, level );

   real_t discr_l2_err = std::sqrt( error.dotGlobal( error, level ) / npoints );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err );

   if ( printTiming )
   {
      walberla::WcTimingTree tt = timingTree->getReduced();
      WALBERLA_LOG_INFO_ON_ROOT( tt );
   }
}

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   test( hyteg::prependHyTeGMeshDir( "quad_8el.msh" ), 4, 20 );
   test( hyteg::prependHyTeGMeshDir( "3D/tet_1el.msh" ), 4, 20 );
   test( hyteg::prependHyTeGMeshDir( "3D/regular_octahedron_8el.msh" ), 4, 20 );

   return 0;
}
