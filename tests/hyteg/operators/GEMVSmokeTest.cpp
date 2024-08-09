/*
* Copyright (c) 2017-2024 Nils Kohl.
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
#include "core/math/Random.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/eigen/EigenSparseMatrix.hpp"
#include "hyteg/eigen/EigenVector.hpp"
#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/mesh/MeshInfo.hpp"
#include "hyteg/petsc/PETScManager.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/primitivestorage/loadbalancing/SimpleBalancer.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesConstantOperator.hpp"
#include "hyteg_operators_composites/stokes/P2P1StokesFullOperator.hpp"

#include "mixed_operator/P2P1TaylorHoodStokesOperator.hpp"

using walberla::real_t;
using namespace hyteg;

template < typename Operator_T >
void compareHytegAndEigenGEMV( const std::string&                         testId,
                               const Operator_T&                          op,
                               const std::shared_ptr< PrimitiveStorage >& storage,
                               const uint_t                               level,
                               bool                                       checkApply = false )
{
   const real_t epsilon = real_c( std::is_same< real_t, double >() ? 1e-15 : 2e-7 );

   using SrcFunction_T          = typename Operator_T::srcType;
   using DstFunction_T          = typename Operator_T::dstType;
   using SrcFunctionNumerator_T = typename SrcFunction_T::template FunctionType< idx_t >;
   using DstFunctionNumerator_T = typename DstFunction_T::template FunctionType< idx_t >;

   SrcFunction_T x( "x", storage, level, level );
   DstFunction_T y( "y", storage, level, level );

   auto alpha = real_c( 0.123 );
   auto beta  = real_c( 0.456 );

   auto rand = []( const Point3D& ) { return walberla::math::realRandom(); };

   x.interpolate( rand, level, All );
   y.interpolate( rand, level, All );

   ///////////////////////
   /// Eigen reference ///
   ///////////////////////

   SrcFunctionNumerator_T xNum( "xNum", storage, level, level );
   DstFunctionNumerator_T yNum( "yNum", storage, level, level );

   DstFunction_T eigenYResult( "eigenYResult", storage, level, level );

   xNum.enumerate( level );
   yNum.enumerate( level );

   Eigen::SparseMatrix< real_t, Eigen::RowMajor > eigenA = createEigenSparseMatrixFromOperator( op, level, xNum, yNum );
   VectorXr                                       eigenX = createEigenVectorFromFunction( x, xNum, level );
   VectorXr                                       eigenY = createEigenVectorFromFunction( y, yNum, level );

   if ( checkApply )
   {
      eigenY = eigenA * eigenX;
   }
   else
   {
      eigenY = alpha * eigenA * eigenX + beta * eigenY;
   }

   assignFunctionFromEigenVector( eigenY, eigenYResult, yNum, level );

   //////////////////
   /// HyTeG GEMV ///
   //////////////////

   if ( checkApply )
   {
      op.apply( x, y, level, All, Replace );
   }
   else
   {
      op.gemv( alpha, x, beta, y, level, All );
   }

   ///////////////////
   /// Error check ///
   ///////////////////

   DstFunction_T error( "error", storage, level, level );

   error.assign( { 1.0, -1.0 }, { y, eigenYResult }, level, All );

   auto errorMaxMagnitude = error.getMaxMagnitude( level );

   auto testString = testId + " " + ( checkApply ? "apply" : "gemv" );
   WALBERLA_LOG_INFO_ON_ROOT( "error inf: " << errorMaxMagnitude << "    | " << testString );
   WALBERLA_CHECK_LESS( errorMaxMagnitude, epsilon );
}

int main( int argc, char* argv[] )
{
   // General setup stuff
   walberla::MPIManager::instance()->initializeMPI( &argc, &argv );
   walberla::MPIManager::instance()->useWorldComm();

   // 3D
   {
      auto                  meshInfo = MeshInfo::fromGmshFile( prependHyTeGMeshDir( "3D/cube_6el.msh" ) );
      SetupPrimitiveStorage setupStorage( meshInfo, walberla::uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      loadbalancing::roundRobin( setupStorage );
      std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

      const uint_t level = 4;

      {
         P1ElementwiseLaplaceOperator op( storage, level, level );
         compareHytegAndEigenGEMV( "P1ElementwiseLaplaceOperator", op, storage, 4 );
         compareHytegAndEigenGEMV( "P1ElementwiseLaplaceOperator", op, storage, 4, true );
      }
   }

   return 0;
}
