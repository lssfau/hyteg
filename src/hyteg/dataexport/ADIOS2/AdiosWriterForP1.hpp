/*
 * Copyright (c) 2023 Marcus Mohr.
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

#include <adios2.h>

#include "hyteg/dataexport/ADIOS2/AdiosWriter.hpp"
#include "hyteg/dataexport/FEFunctionRegistry.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

class AdiosWriter;

class AdiosWriterForP1
{
 public:
   /// \param adios
   /// \param dir             Directory where the files are stored
   /// \param filename        Basename of the vtk files
   /// \param storage         PrimitiveStorage containing the functions
   /// \param writeFrequency  Specifies the frequency of the VTK output see write()
   AdiosWriterForP1( std::string&                               filePath,
                     std::string&                               fileBaseName,
                     std::string&                               configFile,
                     uint_t                                     level,
                     const std::shared_ptr< PrimitiveStorage >& storage );

   /// The destructor takes care of closing the ADIOS2 output file
   ~AdiosWriterForP1()
   {
      WALBERLA_LOG_INFO_ON_ROOT( "D'tor of AdiosWriterForP1 called" );
      // if ( firstWriteCompleted_ )
      {
         engine_.Close();
      }
   }

   /// Trigger exporting data of FE functions of type P1
   ///
   /// \note The caller needs to make sure that the functions have been synced before
   ///       invoking this method!
   void write( const FEFunctionRegistry& registry, uint_t timestep );

 private:
   /// Store the mesh on which our functions live in the output file
   void writeMesh( const std::vector< std::string >& p1FunctionList ) const;

   /// central ADIOS2 interface objects (mutable because write() is marked as const !)
   mutable adios2::ADIOS  adios_;
   mutable adios2::IO     io_;
   mutable adios2::Engine engine_;

   /// storage associated with the P1 type functions
   std::shared_ptr< PrimitiveStorage > storage_{ nullptr };

   /// mesh level on which this writer performs data exports
   uint_t level_;

   /// name of the output file
   std::string fileName_;

   /// need to keep track of chronology (better use TIME later on?)
   mutable bool firstWriteCompleted_{ false };

   /// associate ADIOS variables with P1 type functions in the registry
   ///
   /// \note function will abort, if it encounters and already defined variable
   template < typename value_t >
   void defineVariables( const FEFunctionRegistry& registry ) const;

   template < typename value_t >
   void scheduleScalarFunctionForExport( const P1Function< value_t >& func ) const;

   template < typename value_t >
   void scheduleVectorFunctionForExport( const P1VectorFunction< value_t >& func ) const;
};

} // namespace hyteg
