/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include "SQL.hpp"

namespace hyteg {

void FixedSizeSQLDB::writeRowOnRoot()
{
   for ( auto key : variableKeys_ )
   {
      auto count = variableEntriesInteger_.count( key );
      count += variableEntriesString_.count( key );
      count += variableEntriesDouble_.count( key );

      WALBERLA_CHECK_EQUAL( count, 1, "Variable key \"" << key << "\" not set for this row." );
   }

   std::map< std::string, std::string > stringEntries;
   std::map< std::string, double >      doubleEntries;
   std::map< std::string, int64_t >     integerEntries;

   stringEntries.insert( constantEntriesString_.begin(), constantEntriesString_.end() );
   stringEntries.insert( variableEntriesString_.begin(), variableEntriesString_.end() );

   doubleEntries.insert( constantEntriesDouble_.begin(), constantEntriesDouble_.end() );
   doubleEntries.insert( variableEntriesDouble_.begin(), variableEntriesDouble_.end() );

   integerEntries.insert( constantEntriesInteger_.begin(), constantEntriesInteger_.end() );
   integerEntries.insert( variableEntriesInteger_.begin(), variableEntriesInteger_.end() );

   WALBERLA_ROOT_SECTION()
   {
      db_.storeRun( integerEntries, stringEntries, doubleEntries );
   }

   columnsFixed_ = true;
   variableEntriesInteger_.clear();
   variableEntriesString_.clear();
   variableEntriesDouble_.clear();
}

} // namespace hyteg