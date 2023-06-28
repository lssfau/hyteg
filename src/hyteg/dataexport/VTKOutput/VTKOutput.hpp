/*
 * Copyright (c) 2017-2023 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include <map>
#include <string>
#include <utility>
#include <vector>

#include "core/DataTypes.h"

#include "hyteg/dataexport/FEFunctionRegistry.hpp"
#include "hyteg/dataexport/FEFunctionWriter.hpp"

// our friends and helpers

// clang off
// ordering matters here, otherwise we need to add forward declarations
#include "hyteg/dataexport/VTKOutput/VTKHelpers.hpp"
// clang on

#include "hyteg/dataexport/VTKOutput/VTKEdgeDoFWriter.hpp"
#include "hyteg/dataexport/VTKOutput/VTKMeshWriter.hpp"
#include "hyteg/dataexport/VTKOutput/VTKN1E1Writer.hpp"
#include "hyteg/dataexport/VTKOutput/VTKP1DGEWriter.hpp"
#include "hyteg/dataexport/VTKOutput/VTKP1Writer.hpp"
#include "hyteg/dataexport/VTKOutput/VTKP2Writer.hpp"

// from walblera
#include "vtk/Base64Writer.h"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint64_t;
using walberla::uint_c;
using walberla::uint_t;

class PrimitiveStorage;

class VTKOutput : public FEFunctionWriter
{
 public:
   ///
   /// \param dir             Directory where the files are stored
   /// \param filename        Basename of the vtk files
   /// \param storage         PrimitiveStorage containing the functions
   /// \param writeFrequency  Specifies the frequency of the VTK output see write()
   VTKOutput( std::string                                dir,
              std::string                                filename,
              const std::shared_ptr< PrimitiveStorage >& storage,
              const uint_t&                              writeFrequency = 1 );

   /// Add an FE Function to became part of the next dataexport phase
   template < template < typename > class func_t, typename value_t >
   inline void add( const func_t< value_t >& function )
   {
      feFunctionRegistry_.add< func_t, value_t >( function );
   }

   /// Writes the VTK output only if writeFrequency > 0 and timestep % writeFrequency == 0.
   /// Therefore always writes output if timestep is 0.
   /// Appends the time step to the filename.
   /// Note: files will be overwritten if called twice with the same time step!
   void write( const uint_t& level, const uint_t& timestep = 0 ) const override final;

   /// Set parameter specified by string key to value specified by string value
   ///
   /// The only key currently supported by VTKOutput is "vtkDataFormat" with the two possible values
   /// - ASCII
   /// - BINARY
   void setParameter( const std::string& key, const std::string& value ) override final
   {
      if ( key != "vtkDataFormat" )
      {
         WALBERLA_ABORT( "VTKOutput::setParameter() does not support key = '" << key << "'!" );
      }
      else if ( value == "ASCII" )
      {
         setVTKDataFormat( vtk::DataFormat::ASCII );
      }
      else if ( value == "BINARY" )
      {
         setVTKDataFormat( vtk::DataFormat::BINARY );
      }
      else
      {
         WALBERLA_ABORT( "VTKOutput::setParameter() key vtkDataFormat = '" << value << "' is not supported!" );
      }
   };

   void setVTKDataFormat( vtk::DataFormat vtkDataFormat ) { vtkDataFormat_ = vtkDataFormat; }

 private:
   /// Wrapper class that handles writing data in ASCII or binary format.
   ///
   /// \tparam DTypeInVTK data type that the input data is converted to before writing it to the VTK file
   template < typename DTypeInVTK >
   class VTKStreamWriter
   {
    public:
      explicit VTKStreamWriter( vtk::DataFormat vtkDataFormat )
      : vtkDataFormat_( vtkDataFormat )
      {
         if ( vtkDataFormat_ == vtk::DataFormat::ASCII )
         {
            outputAscii_ << std::scientific;
         }
      }

      template < typename T >
      VTKStreamWriter& operator<<( const T& data )
      {
         if ( vtkDataFormat_ == vtk::DataFormat::ASCII )
         {
            outputAscii_ << static_cast< DTypeInVTK >( data ) << "\n";
         }
         else if ( vtkDataFormat_ == vtk::DataFormat::BINARY )
         {
            outputBase64_ << static_cast< DTypeInVTK >( data );
         }

         return *this;
      }

      void toStream( std::ostream& os )
      {
         if ( vtkDataFormat_ == vtk::DataFormat::ASCII )
         {
            os << outputAscii_.str();
            // reset string stream
            // outputAscii_.str( std::string() );
            // outputAscii_.clear();
         }
         else if ( vtkDataFormat_ == vtk::DataFormat::BINARY )
         {
            outputBase64_.toStream( os );
            // Base64Writer::toStream() already reset the object
            // so nothing left to do for us here
         }
      }

    private:
      vtk::DataFormat             vtkDataFormat_;
      std::ostringstream          outputAscii_;
      walberla::vtk::Base64Writer outputBase64_;
   };

   static const std::map< vtk::DoFType, std::string > DoFTypeToString_;

   void   writeDoFByType( std::ostream& output, const uint_t& level, const vtk::DoFType& dofType ) const;
   uint_t getNumRegisteredFunctions( const vtk::DoFType& dofType ) const;

   std::string fileNameExtension( const vtk::DoFType& dofType, const uint_t& level, const uint_t& timestep ) const;

   void syncAllFunctions( const uint_t& level ) const;

   /// Writes only macro-faces.
   void set2D() { write2D_ = true; }

   /// Writes only macro-cells.
   void set3D() { write2D_ = false; }

   std::string dir_;
   std::string filename_;

   const std::string defaultFMT_ = R"(format="ascii")";

   uint_t writeFrequency_;

   bool write2D_;

   FEFunctionRegistry feFunctionRegistry_;

   std::shared_ptr< PrimitiveStorage > storage_;

   vtk::DataFormat vtkDataFormat_;

   // all writers currently need to be our friends
   friend class VTKFaceDoFWriter;
   friend class VTKEdgeDoFWriter;
   friend class VTKMeshWriter;
   friend class VTKP1Writer;
   friend class VTKP2Writer;
   friend class VTKDGWriter;
   friend class VTKP1DGEWriter;
   friend class VTKN1E1Writer;
};

} // namespace hyteg
