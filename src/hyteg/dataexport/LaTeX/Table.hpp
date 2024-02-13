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

namespace hyteg {
namespace latex {

/// \brief Export tabular material to whitespace separated text files.
///
/// Tables in this format can be directly loaded with pgfplots.
/// This class is templated by the number of columns.
///
/// # Example
///
/// ```cpp
/// Table<2> table( { "Level", "L2error" } );
///
/// table.addElement( 0, 0, 0 );
/// table.addElement( 0, 1, 1e-1 );
///
/// table.addElement( 1, 0, 1 );
/// table.addElement( 1, 1, 2.5e-2 );
///
/// table.write("out", "table");
/// ```
///
/// ```latex
/// \usepackage{booktabs}
/// \usepackage{pgfplotstable}
/// \pgfplotstableset{ every head row/.style = { before row = \toprule
///                                            , after row  = \midrule
///                                            }
///                  , every last row/.style = { after row  = \bottomrule }
///                  }
///
/// \begin{table}
///   \pgfplotstabletypeset[ columns/L2error/.style = { column name = {$\|err\|_{L^2}$} }
///                        ]{out/table.dat}
/// \end{table}
/// ```
template < std::size_t N >
class Table
{
 private:
   std::vector< std::array< std::string, N > > rows_;
   walberla::uint_t                            currentRow_;
   std::stringstream                           stringStream_;

   /// Auxilliary method for use by pushRow()
   template < typename T, class... Args >
   void addRowElement( size_t rowIdx, T& firstArg, Args... args )
   {
      size_t colIdx = N - 1 - sizeof...( Args );
      addElement( rowIdx, colIdx, firstArg );
      if constexpr ( sizeof...( Args ) > 0 )
      {
         addRowElement( rowIdx, args... );
      }
   }

   /// Auxilliary method for printing a specific rows without increasing the current row counter
   void printSpecificRow( std::ostream& os, const std::array<std::basic_string<char>, N> row ) const
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

   /// Auxilliary method for printing the next row that has not been printed with this function yet
   void printNextRow( std::ostream& os )
   {
      WALBERLA_ASSERT_LESS( currentRow_, rows_.size() );
      const auto& row = rows_[currentRow_++];
      printSpecificRow( os, row );
   }

 public:
   /// \brief Create a new `Table` with the given column `headers`.
   Table( std::array< std::string, N >&& headers )
   : rows_{ headers }, currentRow_{ 0 }
   {}

   /// \brief Inserts an element into this table.
   ///
   /// Additional rows are added to the table if necessary.
   /// `T` must support the `stringstream << T` operator.
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

   /// \brief Write this table in whitespace separated format to `dir/filename.dat`.
   void write( const std::string& dir, const std::string& filename ) const
   {
      WALBERLA_ROOT_SECTION()
      {
         std::string   datFilename( walberla::format( "%s/%s.dat", dir.c_str(), filename.c_str() ) );
         std::ofstream file( datFilename );

         file << *this;
      }
   }

   /// \brief Write this table in whitespace separated format to `dir/filename.dat`.
   void writeUpdate( const std::string& dir, const std::string& filename )
   {
      WALBERLA_ROOT_SECTION()
      {
         std::string   datFilename( walberla::format( "%s/%s.dat", dir.c_str(), filename.c_str() ) );
         auto streamMode = ( currentRow_ == 0 ) ? std::ios_base::trunc : std::ios_base::app;
         std::ofstream file( datFilename, streamMode );

         while ( currentRow_ < rows_.size() )
         {
            printNextRow( file );
         }
      }
   }

   /// Append a new row to the end of the table
   template < class... Args >
   void pushRow( Args... args )
   {
      if constexpr ( N != sizeof...( Args ) )
      {
         WALBERLA_ABORT( "Can only use pushRow() with " << N << " arguments!" );
      }

      rows_.emplace_back();
      addRowElement( rows_.size() - 2, args... );
   }

   template < std::size_t M >
   friend std::ostream& operator<<( std::ostream& os, const Table< M >& table );
};

/// \brief Write `table` in whitespace separated format to `os`.
template < std::size_t N >
std::ostream& operator<<( std::ostream& os, const Table< N >& table )
{
   for ( const auto& row : table.rows_ )
   {
      table.printSpecificRow( os, row );
   }

   return os;
}

} // namespace latex
} // namespace hyteg
