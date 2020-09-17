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

#pragma once

#include "core/DataTypes.h"
#include "core/debug/CheckFunctions.h"

#include "sqlite/SQLite.h"

using walberla::double_c;
using walberla::int64_c;
using walberla::int64_t;
using walberla::real_t;
using walberla::uint_t;

namespace hyteg {

/// \brief Wrapper to conveniently export simulation data to SQL, that prohibits empty cells.
class FixedSizeSQLDB
{
 public:
   /// \brief Creates an SQLite database handle.
   ///
   /// \param dbFile path of the file to write the data to
   /// \param allowOverwritingValues if false, values can only be set once per row, guards a little more against accidentally
   ///                               overwriting values
   FixedSizeSQLDB( std::string dbFile, bool allowOverwritingValues = false )
   : db_( dbFile )
   , columnsFixed_( false )
   , allowOverwritingValues_( allowOverwritingValues )
   {}

   /// \brief Adds a constant entry to the database, meaning, this entry is inserted into all rows of following write calls and
   ///        cannot be changed.
   ///@{
   void setConstantEntry( std::string key, std::string value ) { setConstantEntryTpl( key, value, constantEntriesString_ ); }
   void setConstantEntry( std::string key, double value ) { setConstantEntryTpl( key, value, constantEntriesDouble_ ); }
   void setConstantEntry( std::string key, int64_t value ) { setConstantEntryTpl( key, value, constantEntriesInteger_ ); }

   void setConstantEntry( std::string key, uint_t value ) { setConstantEntry( key, int64_c( value ) ); };
   void setConstantEntry( std::string key, bool value ) { setConstantEntry( key, int64_c( value ) ); };
   ///@}

   /// \brief Adds a variable entry to the database, meaning, this entry has to be set before every write call.
   ///@{
   void setVariableEntry( std::string key, std::string value ) { setVariableEntryTpl( key, value, variableEntriesString_ ); }
   void setVariableEntry( std::string key, double value ) { setVariableEntryTpl( key, value, variableEntriesDouble_ ); }
   void setVariableEntry( std::string key, int64_t value ) { setVariableEntryTpl( key, value, variableEntriesInteger_ ); }

   void setVariableEntry( std::string key, uint_t value ) { setVariableEntry( key, int64_c( value ) ); };
   void setVariableEntry( std::string key, bool value ) { setVariableEntry( key, int64_c( value ) ); };
   ///@}

   /// \brief Stores the currently set row into the database file.
   ///
   /// If there is a column entry missing for this row, this method will abort.
   /// This way it is ensured that there are no empty cells.
   void writeRow();

 private:
   template < typename ValueType >
   inline void setConstantEntryTpl( std::string key, ValueType value, std::map< std::string, ValueType >& m )
   {
      WALBERLA_CHECK( !columnsFixed_, "Cannot insert \"" << key << "\" since columns are fixed after first write." )
      if ( !allowOverwritingValues_ )
      {
         WALBERLA_CHECK_EQUAL( m.count( key ), 0, "Attempt to overwrite \"" << key << "\" dismissed." )
      }
      m[key] = value;
   }

   template < typename ValueType >
   void setVariableEntryTpl( std::string key, ValueType value, std::map< std::string, ValueType >& m )
   {
      const bool hasBeenInsertedAnyRow = std::find( variableKeys_.begin(), variableKeys_.end(), key ) != variableKeys_.end();

      if ( columnsFixed_ )
      {
         WALBERLA_CHECK( hasBeenInsertedAnyRow, "New variable key \"" << key << "\" cannot be inserted since columns are fixed." )
      }
      else
      {
         variableKeys_.push_back( key );
      }

      if ( !allowOverwritingValues_ )
      {
         WALBERLA_CHECK_EQUAL( m.count( key ), 0, "Attempt to overwrite \"" << key << "\" dismissed." )
      }


      m[key] = value;
   }

   walberla::sqlite::SQLiteDB db_;
   bool                       columnsFixed_;
   bool                       allowOverwritingValues_;

   std::map< std::string, std::string > constantEntriesString_;
   std::map< std::string, double >      constantEntriesDouble_;
   std::map< std::string, int64_t >     constantEntriesInteger_;

   std::map< std::string, std::string > variableEntriesString_;
   std::map< std::string, double >      variableEntriesDouble_;
   std::map< std::string, int64_t >     variableEntriesInteger_;

   std::vector< std::string > variableKeys_;
};

} // namespace hyteg
