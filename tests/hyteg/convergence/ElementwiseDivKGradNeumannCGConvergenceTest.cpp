/*
 * Copyright (c) 2017-2021 Dominik Thoennes, Nils Kohl.
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSe.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "core/Environment.h"
#include "core/config/Config.h"
#include "core/logging/Logging.h"
#include "core/math/Constants.h"
#include "core/timing/Timer.h"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/forms/form_hyteg_generated/p1/p1_div_k_grad_affine_q3.hpp"
#include "hyteg/forms/form_hyteg_generated/p2/p2_div_k_grad_affine_q4.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/Visualization.hpp"
#include "hyteg/solvers/CGSolver.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

namespace hyteg {

template < typename ElementwiseOperator, typename MassOperator, typename FunctionType, typename DivKGradForm >
void ElementwiseDivKGradCGTest( const uint_t dim, const uint_t level, const real_t targetError )
{
   auto meshInfo = MeshInfo::emptyMeshInfo();
   if ( dim == 2 )
   {
      Point2D n( { 1, 1 } );
      meshInfo = MeshInfo::meshRectangle( -n, n, MeshInfo::CRISS, 2, 2 );
   }
   else
   {
      Point3D n( { 1, 1, 1 } );
      meshInfo = MeshInfo::meshCuboid( -n, n, 1, 1, 1 );
   }
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 2, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "ElementwiseDivKGradNeumannCGConvergenceTest_domain" );

   auto exact = []( const hyteg::Point3D& x ) { return cos( pi * x[0] ) * cos( pi * x[1] ) * cos( pi * x[2] ); };
   auto k     = []( const hyteg::Point3D& x ) { return tanh( 2 * x[0] - 1.0 ) + 2; };
   auto rhs   = [dim]( const hyteg::Point3D& x ) {
      auto x0 = pi * x[0];
      auto x1 = tanh( 2 * x[0] - 1.0 );
      return pi * ( real_t( dim ) * pi * ( x1 + 2 ) * cos( x0 ) - ( 2 * x1 * x1 - 2 ) * sin( x0 ) ) * cos( pi * x[1] ) *
             cos( pi * x[2] );
   };

   DivKGradForm        form( k, k );
   ElementwiseOperator A( storage, level, level, form );
   MassOperator        M( storage, level, level );

   FunctionType r( "r", storage, level, level );
   FunctionType b( "b", storage, level, level );
   FunctionType u( "u", storage, level, level );
   FunctionType u_exact( "u_exact", storage, level, level );
   FunctionType coeff( "k", storage, level, level );
   FunctionType ones( "c", storage, level, level );
   FunctionType err( "err", storage, level, level );

   // initialize discrete functions
   u_exact.interpolate( exact, level );
   coeff.interpolate( k, level );
   r.interpolate( rhs, level );
   M.apply( r, b, level, All );

   // orthogonally project u_exact and b onto range(A)
   ones.interpolate( 1.0, level );
   auto N = ones.dotGlobal( ones, level );
   u_exact.assign( { 1, -ones.dotGlobal( u_exact, level ) / N }, { u_exact, ones }, level );
   b.assign( { 1, -ones.dotGlobal( b, level ) / N }, { b, ones }, level );

   // solve
   auto solver = CGSolver< ElementwiseOperator >( storage, level, level, 1000, 1e-12 );
   solver.solve( A, u, b, level );
   // orthogonally project u onto range(A)
   u.assign( { 1, -ones.dotGlobal( u, level ) / N }, { u, ones }, level );

   // residual and error
   A.apply( u, r, level, All );
   r.assign( { 1.0, -1.0 }, { r, b }, level );
   const real_t rel_res = std::sqrt( r.dotGlobal( r, level ) / b.dotGlobal( b, level ) );
   err.assign( { 1.0, -1.0 }, { u_exact, u }, level );
   M.apply( err, r, level, All );
   const real_t discrete_l2_err = std::sqrt( r.dotGlobal( err, level ) );

   // output
   hyteg::VTKOutput vtkOutput(
       "../../output", "ElementwiseDivKGradNeumannCGConvergenceTest_" + std::to_string( dim ) + "d", storage );
   vtkOutput.add( u );
   vtkOutput.add( u_exact );
   vtkOutput.add( b );
   vtkOutput.add( r );
   vtkOutput.add( coeff );
   vtkOutput.write( level );

   WALBERLA_LOG_INFO_ON_ROOT( "dim: " << dim << ", level: " << level );
   WALBERLA_LOG_INFO_ON_ROOT( "||e||_discrL2 = " << discrete_l2_err );
   WALBERLA_LOG_INFO_ON_ROOT( "||r||/||b||  = " << rel_res );

   WALBERLA_CHECK_LESS( discrete_l2_err, targetError );
}

void runAllTestsP1()
{
   typedef forms::p1_div_k_grad_affine_q3    FormType;
   typedef P1ElementwiseOperator< FormType > ElementwiseOperator;
   typedef P1ElementwiseMassOperator         MassOperator;
   typedef P1Function< real_t >              FunctionType;

   WALBERLA_LOG_INFO_ON_ROOT( "P1 tests" )

   ElementwiseDivKGradCGTest< ElementwiseOperator, MassOperator, FunctionType, FormType >( 2, 3, 6e-2 );
   ElementwiseDivKGradCGTest< ElementwiseOperator, MassOperator, FunctionType, FormType >( 2, 4, 16e-3 );
   ElementwiseDivKGradCGTest< ElementwiseOperator, MassOperator, FunctionType, FormType >( 2, 5, 4e-3 );

   ElementwiseDivKGradCGTest< ElementwiseOperator, MassOperator, FunctionType, FormType >( 3, 3, 6e-1 );
   ElementwiseDivKGradCGTest< ElementwiseOperator, MassOperator, FunctionType, FormType >( 3, 4, 2e-1 );
   ElementwiseDivKGradCGTest< ElementwiseOperator, MassOperator, FunctionType, FormType >( 3, 5, 8e-2 );
}

void runAllTestsP2()
{
   typedef forms::p2_div_k_grad_affine_q4    FormType;
   typedef P2ElementwiseOperator< FormType > ElementwiseOperator;
   typedef P2ElementwiseMassOperator         MassOperator;
   typedef P2Function< real_t >              FunctionType;

   WALBERLA_LOG_INFO_ON_ROOT( "P2 tests" )

   ElementwiseDivKGradCGTest< ElementwiseOperator, MassOperator, FunctionType, FormType >( 2, 2, 5e-3 );
   ElementwiseDivKGradCGTest< ElementwiseOperator, MassOperator, FunctionType, FormType >( 2, 3, 4e-4 );
   if constexpr ( std::is_same_v< real_t, double > )
   {
      // single precision accuracy seems to be not enough
      ElementwiseDivKGradCGTest< ElementwiseOperator, MassOperator, FunctionType, FormType >( 2, 4, 3e-5 );
   }

   ElementwiseDivKGradCGTest< ElementwiseOperator, MassOperator, FunctionType, FormType >( 3, 2, 2e-1 );
   ElementwiseDivKGradCGTest< ElementwiseOperator, MassOperator, FunctionType, FormType >( 3, 3, 3e-2 );
   ElementwiseDivKGradCGTest< ElementwiseOperator, MassOperator, FunctionType, FormType >( 3, 4, 2e-3 );
}

} // namespace hyteg

int main( int argc, char* argv[] )
{
   walberla::Environment walberlaEnv( argc, argv );
   walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
   walberla::MPIManager::instance()->useWorldComm();

   hyteg::runAllTestsP1();
   hyteg::runAllTestsP2();
}