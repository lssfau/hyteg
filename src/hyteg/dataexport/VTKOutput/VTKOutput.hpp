/*
 * Copyright (c) 2017-2025 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include <cstdint>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "core/DataTypes.h"

#include "hyteg/dataexport/FEFunctionWriter.hpp"
#include "hyteg/functions/FEFunctionRegistry.hpp"

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
#include "hyteg/dataexport/VTKOutput/VTKP2PlusBubbleWriter.hpp"
#include "hyteg/dataexport/VTKOutput/VTKP2Writer.hpp"
#include "hyteg/dataexport/VTKOutput/VTKStreamWriter.hpp"

// from walberla
#include "vtk/Base64Writer.h"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint64_t;
using walberla::uint_c;
using walberla::uint_t;

class PrimitiveStorage;

class VTKOutput : public FEFunctionWriter< VTKOutput >
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

   /// Add an FE Function to become part of the next dataexport phase
   template < template < typename > class func_t, typename value_t >
   inline void add( const func_t< value_t >& function )
   {
      // Allowed types for vtk printing
      static_assert( std::is_same_v< value_t, double > || std::is_same_v< value_t, float > ||
                         std::is_same_v< value_t, int32_t > || std::is_same_v< value_t, int64_t >,
                     "The VTK printer is able to print only functions of the types double, float, int32 and int64." );

      // Index vectors of non-nodal FE functions can not be printed directly.
      static_assert( !( (std::is_same_v< value_t, int32_t > || std::is_same_v< value_t, int64_t >) &&(
                         std::is_same_v< func_t< value_t >, DG1Function< value_t > > ||
                         std::is_same_v< func_t< value_t >, dg::DGFunction< value_t > > ||
                         std::is_same_v< func_t< value_t >, dg::DGVectorFunction< value_t > > ||
                         std::is_same_v< func_t< value_t >, n1e1::N1E1VectorFunction< value_t > > ||
                         std::is_same_v< func_t< value_t >, EGFunction< value_t > > ||
                         std::is_same_v< func_t< value_t >, EGP0StokesFunction< value_t > >) ),
                     "You requested to export an integer-valued non-nodal finite element *function*.\n"
                     "Most likely, this is not what you want to do. Presumably, the intent is to print\n"
                     "an index *vector* corresponding to a non-nodal finite element discretization. To\n"
                     "do so, add the degrees of freedoms directly to the `VTKOutput`. For example, use\n"
                     "`VTKOutput::add(*n1e1VectorFunction.getDoFs())`.\n"
                     "Nodal finite element discretizations enjoy the property that the coefficient\n"
                     "vector is exactly the evaluation at the nodes. The VTK printer therefore exports\n"
                     "the coefficient vector directly. On the other hand, functions of non-nodal\n"
                     "discretization must be evaluated first by multiplying the coefficients with the\n"
                     "basis functions. This makes no sense for index vectors." );

      feFunctionRegistry_.add< func_t, value_t >( function );
   }

   /// Remove an FE Function so that it is no longer included in the next dataexport phase
   template < template < typename > class func_t, typename value_t >
   inline void remove( const func_t< value_t >& function )
   {
      feFunctionRegistry_.remove( function );
   }

   /// Writes the VTK output only if writeFrequency > 0 and timestep % writeFrequency == 0.
   /// Therefore always writes output if timestep is 0.
   /// Appends the time step to the filename.
   /// Note: files will be overwritten if called twice with the same time step!
   void write( const uint_t level, const uint_t timestep = 0 );

   /// Set parameter specified by string key to value specified by string value
   ///
   /// The only key currently supported by VTKOutput is "vtkDataFormat" with the two possible values
   /// - ASCII
   /// - BINARY
   void setParameter( const std::string& key, const std::string& value )
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
   static const std::map< vtk::DoFType, std::string > DoFTypeToString_;

   void   writeDoFByType( std::ostream& output, const uint_t& level, const vtk::DoFType& dofType ) const;
   uint_t getNumRegisteredFunctions( const vtk::DoFType& dofType ) const;

   std::string fileNameExtension( const vtk::DoFType& dofType, const uint_t& level, const uint_t& timestep ) const;

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
   friend class VTKP2PlusBubbleWriter;
   friend class VTKDGWriter;
   friend class VTKP1DGEWriter;
   friend class VTKN1E1Writer;
};

} // namespace hyteg
