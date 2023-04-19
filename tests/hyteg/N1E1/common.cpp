/*
* Copyright (c) 2023 Daniel Bauer.
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

#include "common.hpp"

#include "core/debug/CheckFunctions.h"
#include "core/logging/Logging.h"

namespace hyteg {
namespace n1e1 {

/// Test that the L2 convergence factor is ≤ 1/4 (+1/40).
void L2ConvergenceTest( const uint_t                                                                             minLevel,
                        const uint_t                                                                             maxLevel,
                        const System                                                                             system,
                        std::function< real_t( const uint_t level, const System& system, const bool writeVtk ) > test,
                        const bool                                                                               writeVTK )
{
   const real_t l2ConvRate  = 1.0 / 4.0;
   const real_t convRateEps = l2ConvRate * 0.1;
   real_t       err         = test( minLevel, system, writeVTK );

   WALBERLA_LOG_INFO_ON_ROOT( "expected L2 rate: " << l2ConvRate << ", threshold: " << l2ConvRate + convRateEps );
   WALBERLA_LOG_INFO_ON_ROOT( "error level " << minLevel << ": " << std::scientific << err );

   for ( uint_t level = minLevel + 1; level <= maxLevel; level++ )
   {
      const real_t errFiner     = test( level, system, writeVTK );
      const real_t computedRate = errFiner / err;

      WALBERLA_LOG_INFO_ON_ROOT( "error level " << level << ": " << std::scientific << errFiner );
      WALBERLA_LOG_INFO_ON_ROOT( "computed rate level " << level << " / " << level - 1 << ": " << computedRate );

      WALBERLA_CHECK_LESS_EQUAL( computedRate,
                                 l2ConvRate + convRateEps,
                                 "Convergence L2 rate level " << level << " vs level " << level - 1
                                                              << " not sufficiently small (computed: " << computedRate
                                                              << ", estimated + eps: " << l2ConvRate + convRateEps << ")" );
      err = errFiner;
   }
}

} // namespace n1e1
} // namespace hyteg
