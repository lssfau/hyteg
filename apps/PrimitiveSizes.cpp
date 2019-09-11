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
#include "core/DataTypes.h"

#include "hyteg/FunctionProperties.hpp"
#include "hyteg/FunctionTraits.hpp"

using walberla::uint_t;

namespace hyteg {

template< typename FunctionTag_T >
void printPrimitiveSizes()
{
  WALBERLA_LOG_INFO_ON_ROOT( "level |            cell |            face |            edge |          vertex |" );
  WALBERLA_LOG_INFO_ON_ROOT( "------+-----------------+-----------------+-----------------+-----------------+" );
  for ( uint_t level = 2; level < 15; level++ )
  {
    const uint_t vertexSize = numberOfInnerDoFs< FunctionTag_T, Vertex >( level );
    const uint_t edgeSize = numberOfInnerDoFs< FunctionTag_T, Edge >( level );
    const uint_t faceSize = numberOfInnerDoFs< FunctionTag_T, Face >( level );
    const uint_t cellSize = numberOfInnerDoFs< FunctionTag_T, Cell >( level );
    WALBERLA_LOG_INFO_ON_ROOT( std::setw(5) << level << " | " << std::setw(15) << cellSize << " | " << std::setw(15) << faceSize << " | " << std::setw(15) << edgeSize << " | " << std::setw(15) << vertexSize << " |" );
  }
}

}

int main( int argc, char* argv[] )
{

  walberla::Environment walberlaEnv( argc, argv );
  walberla::MPIManager::instance()->useWorldComm();

  WALBERLA_LOG_INFO_ON_ROOT( " --- Primitive Sizes (number of INNER DoFs) --- " );
  WALBERLA_LOG_INFO_ON_ROOT( "P1:" );
  WALBERLA_LOG_INFO_ON_ROOT( "" );
  hyteg::printPrimitiveSizes< hyteg::P1FunctionTag >();
  WALBERLA_LOG_INFO_ON_ROOT( "" );
  WALBERLA_LOG_INFO_ON_ROOT( "P2:" );
  WALBERLA_LOG_INFO_ON_ROOT( "" );
  hyteg::printPrimitiveSizes< hyteg::P2FunctionTag >();
}