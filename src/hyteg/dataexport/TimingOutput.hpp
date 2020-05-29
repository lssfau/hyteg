/*
 * Copyright (c) 2017-2019 Nils Kohl.
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

#pragma once

#include "core/timing/TimingTree.h"
#include "core/timing/TimingJSON.h"

namespace hyteg {

/// \brief Convenience function to synchronize and print a TimingTree.
///
/// Notes: - involves multiple allReduces
///        - the TimingTree is copied - the original TimingTree is not modified
inline void printTimingTree( walberla::WcTimingTree timingTree )
{
   timingTree.synchronize();
   auto ttReduced = timingTree.getReduced().getCopyWithRemainder();
   WALBERLA_LOG_INFO_ON_ROOT( ttReduced );
}

/// \brief Convenience function to synchronize and write a TimingTree to a JSON file.
///
/// Notes: - involves multiple allReduces
///        - the TimingTree is copied - the original TimingTree is not modified
inline void writeTimingTreeJSON( walberla::WcTimingTree timingTree, const std::string& file )
{
   timingTree.synchronize();
   auto           ttReduced = timingTree.getReduced().getCopyWithRemainder();
   WALBERLA_ROOT_SECTION()
   {
     nlohmann::json ttJson;
     walberla::timing::to_json( ttJson, ttReduced );
     std::ofstream jsonOutput;
     jsonOutput.open( file );
     jsonOutput << ttJson.dump( 4 );
     jsonOutput.close();
   }
}


} // namespace hyteg