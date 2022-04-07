/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/composites/P1StokesFunction.hpp"
#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/facedofspace_old/FaceDoFFunction.hpp"
#include "hyteg/functions/BlockFunction.hpp"
#include "hyteg/functions/FunctionMultiStore.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
#include "hyteg/p2functionspace/P2Function.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"

// our friends and helpers

// clang off
// ordering matter here, otherwise we need to add forward declarations
#include "hyteg/dataexport/VTKHelpers.hpp"
// clang on

#include "hyteg/dataexport/VTKFaceDoFWriter.hpp"
#include "hyteg/dataexport/VTKEdgeDoFWriter.hpp"
#include "hyteg/dataexport/VTKMeshWriter.hpp"
#include "hyteg/dataexport/VTKP1Writer.hpp"
#include "hyteg/dataexport/VTKP2Writer.hpp"

// from walblera
#include "vtk/Base64Writer.h"

namespace hyteg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint64_t;
using walberla::uint_c;
using walberla::uint_t;

class PrimitiveStorage;

class VTKOutput
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

   void setVTKDataFormat( vtk::DataFormat vtkDataFormat ) { vtkDataFormat_ = vtkDataFormat; }

   template < typename value_t >
   inline void add( const P0Function< value_t >& function )
   {
      dgFunctions_.push_back( *function.getDGFunction() );
   }

   template < typename value_t >
   inline void add( const P1Function< value_t >& function )
   {
      p1Functions_.push_back( function );
   }

   template < typename value_t >
   inline void add( const P2Function< value_t >& function )
   {
      p2Functions_.push_back( function );
   }

   template < typename value_t >
   inline void add( const P1VectorFunction< value_t >& function )
   {
      p1VecFunctions_.push_back( function );
   }

   template < typename value_t >
   inline void add( const P2VectorFunction< value_t >& function )
   {
      p2VecFunctions_.push_back( function );
   }

   template < typename value_t >
   inline void add( const EdgeDoFFunction< value_t >& function )
   {
      edgeDoFFunctions_.push_back( function );
   }

   template < typename value_t >
   inline void add( const FaceDoFFunction_old< value_t >& function )
   {
      faceDoFFunctions_.push_back( function );
   }

   template < typename value_t >
   inline void add( const dg::DGFunction< value_t >& function )
   {
      dgFunctions_.push_back( function );
   }

   template < typename value_t >
   inline void add( const BlockFunction< value_t >& function )
   {
      for ( uint_t k = 0; k < function.getNumberOfBlocks(); k++ )
      {
         add( function[k] );
      }
   }

   template < typename value_t >
   void add( const GenericFunction< value_t >& function );

   void add( const P1StokesFunction< real_t >& function );
   void add( const P2P1TaylorHoodFunction< real_t >& function );

   /// Writes the VTK output only if writeFrequency > 0 and timestep % writeFrequency == 0.
   /// Therefore always writes output if timestep is 0.
   /// Appends the time step to the filename.
   /// Note: files will be overwritten if called twice with the same time step!
   void write( const uint_t& level, const uint_t& timestep = 0 ) const;

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
         }
         else if ( vtkDataFormat_ == vtk::DataFormat::BINARY )
         {
            outputBase64_.toStream( os );
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

   template < template < typename > class WrapperFunc, typename func_t >
   void unwrapAndAdd( const WrapperFunc< func_t >& function )
   {
      add( function.unwrap() );
   }

   template < typename WrapperFunc, typename value_t >
   bool tryUnwrapAndAdd( const GenericFunction< value_t >& function )
   {
      bool               success = false;
      const WrapperFunc* aux     = dynamic_cast< const WrapperFunc* >( &function );
      if ( aux != nullptr )
      {
         add( aux->unwrap() );
         success = true;
      }
      return success;
   }

   std::string dir_;
   std::string filename_;

   const std::string defaultFMT_ = "format=\"ascii\"";

   uint_t writeFrequency_;

   bool write2D_;

   FunctionMultiStore< P1Function > p1Functions_;
   FunctionMultiStore< P2Function > p2Functions_;

   FunctionMultiStore< P1VectorFunction > p1VecFunctions_;
   FunctionMultiStore< P2VectorFunction > p2VecFunctions_;

   FunctionMultiStore< EdgeDoFFunction > edgeDoFFunctions_;
   FunctionMultiStore< FaceDoFFunction_old > faceDoFFunctions_;

   FunctionMultiStore< dg::DGFunction > dgFunctions_;

   std::shared_ptr< PrimitiveStorage > storage_;

   vtk::DataFormat vtkDataFormat_;

   // all writers currently need to be our friends
   friend class VTKFaceDoFWriter;
   friend class VTKEdgeDoFWriter;
   friend class VTKMeshWriter;
   friend class VTKP1Writer;
   friend class VTKP2Writer;
   friend class VTKDGWriter;
};

} // namespace hyteg
