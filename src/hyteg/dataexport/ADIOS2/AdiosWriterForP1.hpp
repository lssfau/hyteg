/*
 * Copyright (c) 2023-2024 Marcus Mohr.
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
#include "hyteg/functions/FEFunctionRegistry.hpp"

namespace hyteg {

using walberla::real_t;
using walberla::uint_t;

using adiosHelpers::adiostype_t;

class AdiosWriter;

class AdiosWriterForP1
{
 public:
   /// \param adios          top-level ADIOS2 interface object
   /// \param filePath       Path to directory where the files are stored
   /// \param fileBaseName   Basename of the vtk files
   /// \param engineType     for file I/O should by a BP format like "BP4"
   /// \param storage        PrimitiveStorage associated with functions to export
   /// \param level          fixed refinement level associated with the writer object
   AdiosWriterForP1( adios2::ADIOS&                             adios,
                     const std::string&                         filePath,
                     const std::string&                         fileBaseName,
                     const std::string&                         engineType,
                     uint_t                                     level,
                     const std::shared_ptr< PrimitiveStorage >& storage );

   /// The destructor takes care of closing the ADIOS2 output file
   ~AdiosWriterForP1()
   {
      if ( firstWriteCompleted_ )
      {
         engine_.Close();
      }
   }

   /// Trigger exporting data of FE functions of type P1
   ///
   /// \note The caller needs to make sure that the functions have been synced before
   ///       invoking this method!
   void write( const FEFunctionRegistry&            registry,
               uint_t                               timestep,
               adios2::Params&                      userProvidedParameters,
               std::map< std::string, adiostype_t > additionalAttributes );

 private:
   /// Store the mesh on which our functions live in the output file
   void writeMesh( const std::vector< std::string >& p1FunctionList );

   /// central ADIOS2 interface objects specific to this writer object
   adios2::IO     io_;
   adios2::Engine engine_;

   /// storage associated with the P1 type functions
   std::shared_ptr< PrimitiveStorage > storage_{ nullptr };

   /// mesh level on which this writer performs data exports
   uint_t level_;

   /// name of the output file
   std::string fileName_;

   /// need to keep track of chronology
   bool firstWriteCompleted_{ false };

   /// associate ADIOS variables with P1 type functions in the registry
   ///
   /// \note function will abort, if it encounters an already defined variable
   template < typename value_t >
   void defineVariables( const FEFunctionRegistry& registry );

   /// Copy function data into an ADIOS span to schedule it for export
   template < typename value_t >
   void scheduleScalarFunctionForExport( const P1Function< value_t >& func );

   /// Copy function data into an ADIOS span to schedule it for export
   template < typename value_t >
   void scheduleVectorFunctionForExport( const P1VectorFunction< value_t >& func );
};

} // namespace hyteg
