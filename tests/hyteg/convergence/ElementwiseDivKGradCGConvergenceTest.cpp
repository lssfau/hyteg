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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
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

#include "hyteg_operators/operators/div_k_grad/P1ElementwiseDivKGrad.hpp"
#include "hyteg_operators/operators/div_k_grad/P2ElementwiseDivKGrad.hpp"
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

#include "constant_stencil_operator/P2ConstantOperator.hpp"

using walberla::real_t;
using walberla::uint_c;
using walberla::uint_t;
using walberla::math::pi;

namespace hyteg {

template < typename Operator >
Operator operatorFactory( const std::shared_ptr< PrimitiveStorage >&       storage,
                          size_t                                           level,
                          std::function< real_t( const hyteg::Point3D& ) > k );

template <>
P1ElementwiseOperator< forms::p1_div_k_grad_affine_q3 > operatorFactory( const std::shared_ptr< PrimitiveStorage >&       storage,
                                                                         size_t                                           level,
                                                                         std::function< real_t( const hyteg::Point3D& ) > k )
{
   forms::p1_div_k_grad_affine_q3 form( k, k );
   return P1ElementwiseOperator( storage, level, level, form );
}

template <>
P2ElementwiseOperator< forms::p2_div_k_grad_affine_q4 > operatorFactory( const std::shared_ptr< PrimitiveStorage >&       storage,
                                                                         size_t                                           level,
                                                                         std::function< real_t( const hyteg::Point3D& ) > k )
{
   forms::p2_div_k_grad_affine_q4 form( k, k );
   return P2ElementwiseOperator( storage, level, level, form );
}

template <>
operatorgeneration::P1ElementwiseDivKGrad operatorFactory( const std::shared_ptr< PrimitiveStorage >&       storage,
                                                           size_t                                           level,
                                                           std::function< real_t( const hyteg::Point3D& ) > k )
{
   P1Function< real_t > kFunction( "k", storage, level, level );
   kFunction.interpolate( k, level );
   return operatorgeneration::P1ElementwiseDivKGrad( storage, level, level, kFunction );
}

template <>
operatorgeneration::P2ElementwiseDivKGrad operatorFactory( const std::shared_ptr< PrimitiveStorage >&       storage,
                                                           size_t                                           level,
                                                           std::function< real_t( const hyteg::Point3D& ) > k )
{
   P2Function< real_t > kFunction( "k", storage, level, level );
   kFunction.interpolate( k, level );
   return operatorgeneration::P2ElementwiseDivKGrad( storage, level, level, kFunction );
}

template < typename ElementwiseOperator, typename MassOperator, typename FunctionType >
void ElementwiseDivKGradCGTest( const std::string& meshFile, const uint_t level, const real_t targetError )
{
   const auto            meshInfo = MeshInfo::fromGmshFile( meshFile );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   auto storage = std::make_shared< PrimitiveStorage >( setupStorage );
   writeDomainPartitioningVTK( storage, "../../output", "ElementwiseDivKGradCGConvergenceTest_domain" );

   // alpha->inf => k->step
   real_t alpha = 2;
   // increase phi to decrease smoothness of u
   real_t phi = 2;

   std::function< real_t( const hyteg::Point3D& ) > exact = [phi]( const hyteg::Point3D& x ) {
      return sin( phi * pi * x[0] ) * sinh( pi * x[1] );
   };
   std::function< real_t( const hyteg::Point3D& ) > boundary = exact;
   std::function< real_t( const hyteg::Point3D& ) > k        = [alpha]( const hyteg::Point3D& x ) {
      return tanh( alpha * ( x[0] - 0.5 ) ) + 2;
   };
   std::function< real_t( const hyteg::Point3D& ) > rhs = [phi, alpha]( const hyteg::Point3D& x ) {
      real_t t0 = tanh( alpha * ( x[0] - 0.5 ) );
      real_t t1 = phi * alpha * ( t0 * t0 - 1 ) * cos( phi * pi * x[0] );
      real_t t2 = ( phi * phi - 1 ) * pi * ( t0 + 2 ) * sin( phi * pi * x[0] );
      return pi * ( t1 + t2 ) * sinh( pi * x[1] );
   };

   ElementwiseOperator L = operatorFactory< ElementwiseOperator >( storage, level, k );
   MassOperator        M( storage, level, level );

   FunctionType r( "r", storage, level, level );
   FunctionType f( "f", storage, level, level );
   FunctionType u( "u", storage, level, level );
   FunctionType u_exact( "u_exact", storage, level, level );
   FunctionType err( "err", storage, level, level );
   FunctionType coeff( "k", storage, level, level );

   coeff.interpolate( k, level );

   r.interpolate( rhs, level );
   M.apply( r, f, level, All );

   u.interpolate( boundary, level, hyteg::DirichletBoundary );
   u_exact.interpolate( exact, level );

   auto solver = CGSolver< ElementwiseOperator >( storage, level, level, 1000, 1e-12 );
   solver.solve( L, u, f, level );

   err.assign( { 1.0, -1.0 }, { u, u_exact }, level );

   const auto   npoints      = real_c( numberOfGlobalDoFs< P2FunctionTag >( *storage, level ) );
   const real_t discr_l2_err = std::sqrt( err.dotGlobal( err, level ) / npoints );

   hyteg::VTKOutput vtkOutput( "../../output", "ElementwiseDivKGradCGConvergenceTest", storage );
   vtkOutput.add( u );
   vtkOutput.add( u_exact );
   vtkOutput.add( f );
   vtkOutput.add( r );
   vtkOutput.add( err );
   vtkOutput.add( coeff );
   vtkOutput.write( level );

   WALBERLA_LOG_INFO_ON_ROOT( "discrete L2 error = " << discr_l2_err << " (level " << level << ", tol: " << targetError
                                                     << ", mesh: " << meshFile << ")" );
   WALBERLA_CHECK_LESS( discr_l2_err, targetError );
}

void runAllTestsP1()
{
   using MassOperator = P1ElementwiseMassOperator;
   using FunctionType = P1Function< real_t >;

   WALBERLA_LOG_INFO_ON_ROOT( "P1 tests" )

   {
      using Operator = P1ElementwiseOperator< forms::p1_div_k_grad_affine_q3 >;

      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/quad_4el.msh", 3, 4.0e-2 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/quad_4el.msh", 4, 1.0e-2 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/quad_4el.msh", 5, 2.5e-3 );

      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/3D/tet_1el.msh", 3, 1.1e-3 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/3D/tet_1el.msh", 4, 3.5e-4 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/3D/tet_1el.msh", 5, 1.0e-4 );
   }
   {
      using Operator = operatorgeneration::P1ElementwiseDivKGrad;

      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/quad_4el.msh", 3, 4.0e-2 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/quad_4el.msh", 4, 1.0e-2 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/quad_4el.msh", 5, 2.5e-3 );

      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/3D/tet_1el.msh", 3, 1.1e-3 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/3D/tet_1el.msh", 4, 3.5e-4 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/3D/tet_1el.msh", 5, 1.0e-4 );
   }
}

void runAllTestsP2()
{
   using MassOperator = P2ElementwiseMassOperator;
   using FunctionType = P2Function< real_t >;

   WALBERLA_LOG_INFO_ON_ROOT( "P2 tests" )

   {
      using Operator = P2ElementwiseOperator< forms::p2_div_k_grad_affine_q4 >;

      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/quad_4el.msh", 2, 2e-2 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/quad_4el.msh", 3, 2e-3 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/quad_4el.msh", 4, 2e-4 );

      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/3D/tet_1el.msh", 2, 3e-3 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/3D/tet_1el.msh", 3, 3e-4 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/3D/tet_1el.msh", 4, 3e-5 );
   }
   {
      using Operator = operatorgeneration::P2ElementwiseDivKGrad;

      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/quad_4el.msh", 2, 2e-2 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/quad_4el.msh", 3, 2e-3 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/quad_4el.msh", 4, 2e-4 );

      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/3D/tet_1el.msh", 2, 3e-3 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/3D/tet_1el.msh", 3, 3e-4 );
      ElementwiseDivKGradCGTest< Operator, MassOperator, FunctionType >( "../../meshes/3D/tet_1el.msh", 4, 3e-5 );
   }
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
