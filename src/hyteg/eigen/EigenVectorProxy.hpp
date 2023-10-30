/*
* Copyright (c) 2023 Nils Kohl.
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

#include "hyteg/eigen/EigenWrapper.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"
#include "hyteg/types/Matrix.hpp"

namespace hyteg {

using walberla::uint_t;

class EigenConstVectorProxy : public VectorProxy
{
 public:
   /// \brief Create a read-only view/proxy onto the passed vector.
   ///
   /// The caller must ensure that `vectorRef` outlives `this`.
   EigenConstVectorProxy( const VectorXr& vectorRef )
   : vectorRef_( vectorRef )
   {}

   /// \brief Sets the passed value in the vector.
   virtual void setValue( uint_t idx, real_t value )
   {
      WALBERLA_ABORT( "Use the non-const version of the proxy to set values." );
   }

   /// \brief Returns the passed value of the vector.
   virtual real_t getValue( uint_t idx ) const { return vectorRef_( static_cast< Eigen::Index >( idx ) ); }

 private:
   const VectorXr& vectorRef_;
};

class EigenVectorProxy : public VectorProxy
{
 public:
   /// \brief Create a mutable view/proxy onto the passed vector.
   ///
   /// The caller must ensure that `vectorRef` outlives `this`.
   EigenVectorProxy( VectorXr& vectorRef )
   : vectorRef_( vectorRef )
   {}

   /// \brief Sets the passed value in the vector.
   virtual void setValue( uint_t idx, real_t value ) { vectorRef_( static_cast< Eigen::Index >( idx ) ) = value; }

   /// \brief Returns the passed value of the vector.
   virtual real_t getValue( uint_t idx ) const { return vectorRef_( static_cast< Eigen::Index >( idx ) ); }

 private:
   VectorXr& vectorRef_;
};

} // namespace hyteg
