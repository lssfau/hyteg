/*
 * Copyright (c) 2017-2019 Dominik Thoennes.
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
// test that the product one^T*M*one with mass matrix M and vector of ones gives area of domain
#include "hyteg/dataexport/VTKOutput.hpp"
#include "hyteg/elementwiseoperators/P2ElementwiseOperator.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/p2functionspace/P2ConstantOperator.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "core/Environment.h"

using walberla::real_t;
using walberla::uint_t;
using namespace hyteg;

typedef enum { P2CONSTANT, P2ELEMENTWISE } opType;

template<typename OperatorType>
void checkArea( std::shared_ptr<PrimitiveStorage> storage, real_t area )
{

  const size_t minLevel = 2;
  const size_t maxLevel = 4;

  OperatorType massOp( storage, minLevel, maxLevel );

  // Some overhead currently required to get the elementwise case correct due to the
  // fact that boundaryTypeToSkipDuringAdditiveCommunication is property of P2Function
  P2Function< real_t > aux( "aux", storage, minLevel, maxLevel, BoundaryCondition::create012BC(), DoFType::None );
  // P2Function< real_t > aux( "aux", storage, minLevel, maxLevel );
  P2Function< real_t > vecOfOnes( "vecOfOnes", storage, minLevel, maxLevel );

  for( uint_t lvl = minLevel; lvl <= maxLevel; ++lvl )
    {
      vecOfOnes.interpolate( real_c(1.0), lvl, All );
      massOp.apply( vecOfOnes, aux, lvl, All );
      real_t measure = vecOfOnes.dotGlobal( aux, lvl );
      WALBERLA_LOG_INFO_ON_ROOT( "measure = " << std::scientific << measure );
      WALBERLA_CHECK_FLOAT_EQUAL( measure, area );
    }
}


int main(int argc, char **argv)
{
  walberla::debug::enterTestMode();

  walberla::mpi::Environment MPIenv(argc, argv);
  walberla::logging::Logging::instance()->setLogLevel( walberla::logging::Logging::PROGRESS );
  walberla::MPIManager::instance()->useWorldComm();

  // Test with rectangle
  WALBERLA_LOG_INFO_ON_ROOT( "Testing with RECTANGLE" );
  MeshInfo meshInfo = MeshInfo::meshRectangle( Point2D( {0.0, -1.0} ), Point2D( {2.0, 3.0} ),
                                               MeshInfo::CRISSCROSS, 1, 2 );
  SetupPrimitiveStorage setupStorage(meshInfo, uint_c(walberla::mpi::MPIManager::instance()->numProcesses()));
  std::shared_ptr<PrimitiveStorage> storage = std::make_shared<PrimitiveStorage>(setupStorage);
  checkArea<P2ConstantMassOperator>( storage, 8.0 );
  checkArea<P2ElementwiseMassOperator>( storage, 8.0 );

  // Test with backward facing step
  WALBERLA_LOG_INFO_ON_ROOT( "Testing with BFS" );
  meshInfo = MeshInfo::fromGmshFile( "../../data/meshes/bfs_12el.msh" );
  SetupPrimitiveStorage setupStorageBFS( meshInfo,
                                         uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
  std::shared_ptr<PrimitiveStorage> storageBFS = std::make_shared<PrimitiveStorage>( setupStorageBFS );
  checkArea<P2ConstantMassOperator>( storageBFS, 1.75 );
  checkArea<P2ElementwiseMassOperator>( storageBFS, 1.75 );

  return EXIT_SUCCESS;
}
