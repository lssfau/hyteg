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

#include <array>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "core/Format.hpp"
#include "core/mpi/MPIManager.h"

template < std::size_t N >
class Table
{
 private:
   std::vector< std::array< std::string, N > > rows_;
   std::stringstream                           stringStream_;

 public:
   Table( std::array< std::string, N >&& headers )
   : rows_{ headers }
   {}

   template < typename T >
   void addElement( const size_t row, const size_t col, const T& elem )
   {
      stringStream_.str( "" );
      stringStream_ << elem;

      if ( row + 1 >= rows_.size() )
      {
         rows_.resize( row + 2 );
      }
      rows_[row + 1][col] = stringStream_.str();
   }

   void write( const std::string& dir, const std::string& filename )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::string   datFilename( walberla::format( "%s/%s.dat", dir.c_str(), filename.c_str() ) );
         std::ofstream file( datFilename );

         file << *this;
      }
   }

   template < std::size_t M >
   friend std::ostream& operator<<( std::ostream& os, const Table< M >& table );
};

template < std::size_t N >
std::ostream& operator<<( std::ostream& os, const Table< N >& table )
{
   for ( const auto& row : table.rows_ )
   {
      for ( auto it = row.begin(); it != row.end(); ++it )
      {
         if ( it != row.begin() )
         {
            os << " ";
         }
         os << *it;
      }
      os << "\n";
   }

   return os;
}
