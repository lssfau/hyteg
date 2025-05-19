/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Marcus Mohr.
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

#include "hyteg/p0functionspace/P0Function.hpp"

#include "core/Environment.h"
#include "core/math/Constants.h"
#include "core/math/Random.h"

#include "hyteg/geometry/AffineMap2D.hpp"
#include "hyteg/geometry/AnnulusMap.hpp"
#include "hyteg/geometry/IcosahedralShellMap.hpp"
#include "hyteg/primitivestorage/SetupPrimitiveStorage.hpp"
#include "hyteg/communication/Syncing.hpp"
#include "hyteg/gridtransferoperators/P1toP0Conversion.hpp"
#include "hyteg/gridtransferoperators/P0toP0AveragedInjection.hpp"

#include "hyteg/dataexport/VTKOutput/VTKOutput.hpp"

using namespace hyteg;

template < typename BlendingMap >
void averagingTestAnnulusSphShell( uint_t nTan, uint_t nRad, uint_t minLevel, uint_t maxLevel, real_t rMin, real_t rMax,
std::vector< real_t > tolLevels )
{
   std::shared_ptr< SetupPrimitiveStorage > setupStorage;
   std::shared_ptr< PrimitiveStorage > storage;

   real_t volFormula;

   if constexpr (std::is_same_v< BlendingMap, AnnulusMap >)
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Running test with AnnulusMap" );

      MeshInfo meshInfo = MeshInfo::meshAnnulus(rMin, rMax, MeshInfo::CRISS, nTan, nRad);
      setupStorage = std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      AnnulusMap::setMap(*setupStorage);
      storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );

      volFormula = (walberla::math::pi / 2.0) * (std::pow(rMax, 4) - std::pow(rMin, 4));
   }
   else if constexpr (std::is_same_v< BlendingMap, IcosahedralShellMap >)
   {
      WALBERLA_LOG_INFO_ON_ROOT( "Running test with IcosahedralShellMap" );

      MeshInfo meshInfo = MeshInfo::meshSphericalShell(nTan, nRad, rMin, rMax);
      setupStorage = std::make_shared< SetupPrimitiveStorage >( meshInfo, uint_c( walberla::mpi::MPIManager::instance()->numProcesses() ) );
      IcosahedralShellMap::setMap(*setupStorage);
      storage = std::make_shared< PrimitiveStorage >( *setupStorage, 1 );

      volFormula = (4.0 * walberla::math::pi / 5.0) * (std::pow(rMax, 5) - std::pow(rMin, 5));
   }
   else
   {
      WALBERLA_ABORT("Called with unsupported blending map");
   }   

   P1Function< real_t > TP1("TP1", storage, minLevel, maxLevel);
   P0Function< real_t > T("T", storage, minLevel, maxLevel);
   P0Function< real_t > cellVol("cellVol", storage, minLevel, maxLevel);

   std::function< real_t(const Point3D&) > TInterp = [&](const Point3D& x)
   {
      real_t r = x.norm();
      return r * r;
   };

   TP1.interpolate(TInterp, maxLevel, All);

   communication::syncFunctionBetweenPrimitives(TP1, maxLevel);

   P0toP0AveragedInjection p0Top0AveragedInjection( hyteg::AveragingType::ARITHMETIC, false );

   // T.averageFromP1(TP1, maxLevel, p0averaging::AVERAGING_METHOD::ARITHMETIC);
   // T.transferToAllLowerLevels(maxLevel, p0averaging::AVERAGING_METHOD::ARITHMETIC, false);

   P1toP0Conversion(TP1, T, maxLevel, hyteg::AveragingType::ARITHMETIC);
   p0Top0AveragedInjection.restrictToAllLowerLevels( T, maxLevel );

   // WALBERLA_LOG_INFO_ON_ROOT("volFormula = " << volFormula);

   for(uint_t level = maxLevel; level > minLevel; level--)
   {
      cellVol.writeElementVolumesToDoFs(level);
      real_t volCalculated = T.dotGlobal(cellVol, level);

      real_t relError = std::abs(volCalculated - volFormula) / volFormula;

      // WALBERLA_LOG_INFO_ON_ROOT("Error at level " << level << " = " << relError );
      WALBERLA_CHECK_LESS( relError, tolLevels[level] );
   }
}

int main( int argc, char** argv )
{
   walberla::debug::enterTestMode();
   walberla::mpi::Environment MPIenv( argc, argv );
   walberla::MPIManager::instance()->useWorldComm();

   averagingTestAnnulusSphShell< AnnulusMap >(8u, 4u, 0u, 5u, 1.22, 2.22, {0.0, 3e-2, 7e-3, 2e-3, 4e-4, 1e-4});
   averagingTestAnnulusSphShell< IcosahedralShellMap >(3u, 2u, 0u, 4u, 1.22, 2.22, {0.0, 4e-2, 8e-3, 2e-3, 3e-4});

   return EXIT_SUCCESS;
}
