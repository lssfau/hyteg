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

#pragma once

#include <fstream>
#include <iomanip>
#include <map>
#include <ostream>
#include <sstream>
#include <string>

#include "core/Format.hpp"
#include "core/mpi/MPIManager.h"

namespace hyteg {

/// \brief Export unstructured data in a LaTeX friendly format.
///
/// # Example
///
/// ```cpp
/// KeyValueStore store;
///
/// store.store( "/example/minLevel", minLevel );
/// store.store( "/example/maxLevel", maxLevel );
///
/// store.writePgfKeys( "out", "example" );
/// ```
///
/// ```latex
/// \usepackage{pgfkeys}
/// \usepackage{siunitx}
///
/// \input{out/example.tex}
/// The coarsest and finest grids are on level \num{\pgfkeysvalueof{/example/minLevel}} and \num{\pgfkeysvalueof{/example/maxLevel}}, respectively.
/// ```
class KeyValueStore
{
 private:
   std::map< std::string, std::string > keyValues_;

   std::size_t       maxKeyLength_ = 0;
   std::stringstream stringStream_;

 public:
   /// \brief Store a single key-value pair.
   ///
   /// `T` must support the `stringstream << T` operator.
   template < typename T >
   void store( std::string key, T value )
   {
      stringStream_.str( "" );
      stringStream_ << value;
      keyValues_[key] = stringStream_.str();

      maxKeyLength_ = std::max( maxKeyLength_, key.length() );
   }

   /// \brief Write a tex-file defining pgfkeys key-value pairs to `dir/filename.tex`.
   void writePgfKeys( const std::string& dir, const std::string& filename )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::string   texFilename( walberla::format( "%s/%s.tex", dir.c_str(), filename.c_str() ) );
         std::ofstream file( texFilename );

         for ( auto it = keyValues_.begin(); it != keyValues_.end(); ++it )
         {
            file << "\\pgfkeyssetvalue{" << it->first << "}{" << it->second << "}\n";
         }
      }
   }

   friend std::ostream& operator<<( std::ostream& os, const KeyValueStore& store );
};

/// \brief Write all key-value pairs in a human-readable format to `os`.
inline std::ostream& operator<<( std::ostream& os, const KeyValueStore& store )
{
   for ( auto it = store.keyValues_.begin(); it != store.keyValues_.end(); ++it )
   {
      os << std::left << std::setw( static_cast< int >( store.maxKeyLength_ ) ) << it->first << " = " << it->second << "\n";
   }
   return os;
}

} // namespace hyteg
