/*
 * Copyright (c) 2020 Andreas Wagner
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
#include "core/math/Constants.h"
#include "core/mpi/MPIManager.h"

#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"

#include "hyteg/elementwiseoperators/P1ElementwiseOperator.hpp"
#include "hyteg/p1functionspace/P1ConstantOperator.hpp"
#include "hyteg/p1functionspace/P1VariableOperator.hpp"

#include "hyteg/solvers/GaussSeidelSmoother.hpp"
#include "hyteg/solvers/WeightedJacobiSmoother.hpp"
#include "hyteg/solvers/SymmetricGaussSeidelSmoother.hpp"
#include "hyteg/solvers/ChebyshevSmoother.hpp"
#include "hyteg/solvers/SORSmoother.hpp"
#include "hyteg/solvers/SymmetricSORSmoother.hpp"

using walberla::real_t;
using walberla::math::pi;

using namespace hyteg;

template < typename SmootherType, typename OperatorType, typename FunctionType  >
bool smootherThrowsException(SmootherType& smoother, OperatorType& op, FunctionType& src, FunctionType& dst, uint_t level){
   try {
      smoother.solve( op, dst, src, level );
   } catch (const std::runtime_error& e) {
      return true;
   }
   return false;
}

int main( int argc, char** argv )
{
   walberla::Environment env( argc, argv );
   walberla::mpi::MPIManager::instance()->useWorldComm();

   const uint_t minLevel = 2;
   const uint_t maxLevel = 3;

   MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( { -1, -1 } ), Point2D( { 1., 1. } ), MeshInfo::CRISSCROSS, 2, 2 );
   SetupPrimitiveStorage setupStorage( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
   setupStorage.setMeshBoundaryFlagsOnBoundary( 1, 0, true );
   std::shared_ptr< PrimitiveStorage > storage = std::make_shared< PrimitiveStorage >( setupStorage );

   P1Function< real_t > p1Src( "p1Src", storage, minLevel, maxLevel );
   P1Function< real_t > p1Dst( "p1Dst", storage, minLevel, maxLevel );

   using P1ElOp = P1ElementwiseLaplaceOperator;
   using P1ConstOp = P1ConstantLaplaceOperator;
   using P1VarOp = P1BlendingLaplaceOperator;

   P1ElOp p1ElOp(storage, minLevel, maxLevel);
   p1ElOp.computeInverseDiagonalOperatorValues();
   P1ConstOp p1ConstOp(storage, minLevel, maxLevel);
   P1VarOp p1VarOp(storage, minLevel, maxLevel);
   p1VarOp.computeInverseDiagonalOperatorValues();

   GaussSeidelSmoother< P1ElOp > gsEl;
   WeightedJacobiSmoother< P1ElOp > wJacEl(storage, minLevel, maxLevel, maxLevel);
   ChebyshevSmoother< P1ElOp > chebEl(storage, minLevel, maxLevel);
   chebEl.setupCoefficients(1, 1);
   WALBERLA_CHECK( !smootherThrowsException(wJacEl, p1ElOp, p1Src, p1Dst, minLevel) );
   WALBERLA_CHECK( !smootherThrowsException(chebEl, p1ElOp, p1Src, p1Dst, minLevel) );
   WALBERLA_CHECK( smootherThrowsException(gsEl, p1ElOp, p1Src, p1Dst, minLevel) );

   GaussSeidelSmoother< P1ConstOp > gsConst;
   SymmetricGaussSeidelSmoother< P1ConstOp > sgsConst;
   SymmetricSORSmoother< P1ConstOp > ssorConst(1.);
   SORSmoother< P1ConstOp > sorConst(1.);
   WALBERLA_CHECK( !smootherThrowsException(gsConst, p1ConstOp, p1Src, p1Dst, minLevel) );
   WALBERLA_CHECK( !smootherThrowsException(sgsConst, p1ConstOp, p1Src, p1Dst, minLevel) );
   WALBERLA_CHECK( !smootherThrowsException(ssorConst, p1ConstOp, p1Src, p1Dst, minLevel) );
   WALBERLA_CHECK( !smootherThrowsException(sorConst, p1ConstOp, p1Src, p1Dst, minLevel) );

   ChebyshevSmoother< P1VarOp > chebVar(storage, minLevel, maxLevel);
   chebVar.setupCoefficients(1, 1);
   SymmetricGaussSeidelSmoother< P1VarOp > sgsVar;
   SymmetricSORSmoother< P1VarOp > ssorVar(1.);
   SORSmoother< P1VarOp > sorVar(1.);
   WeightedJacobiSmoother< P1VarOp > wJacVar(storage, minLevel, maxLevel, maxLevel);
   WALBERLA_CHECK( !smootherThrowsException(chebVar, p1VarOp, p1Src, p1Dst, minLevel) );
   WALBERLA_CHECK( smootherThrowsException(sgsVar, p1VarOp, p1Src, p1Dst, minLevel) );
   WALBERLA_CHECK( smootherThrowsException(ssorVar, p1VarOp, p1Src, p1Dst, minLevel) );
   WALBERLA_CHECK( !smootherThrowsException(sorVar, p1VarOp, p1Src, p1Dst, minLevel) );
   WALBERLA_CHECK( !smootherThrowsException(wJacVar, p1VarOp, p1Src, p1Dst, minLevel) );
}
