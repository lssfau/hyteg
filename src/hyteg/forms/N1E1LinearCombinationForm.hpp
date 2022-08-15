/*
 * Copyright (c) 2022 Daniel Bauer.
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

#include <vector>

#include "core/debug/CheckFunctions.h"

#include "hyteg/forms/N1E1Form.hpp"

namespace hyteg {
namespace n1e1 {

/// \brief N1E1 form that performs a linear combination of other N1E1 forms.
///
/// The linear combinations are formed according to the passed scalars and forms:
///    result = scalar[0] * form[0]->integrate() + scalar[1] * form[1]->integrate() + ...
/// The number of passed scalars and forms must be identical.
///
class N1E1LinearCombinationForm : public N1E1Form
{
 public:
   N1E1LinearCombinationForm() = default;

   N1E1LinearCombinationForm( const std::vector< real_t >& scalars, const std::vector< N1E1Form* >& forms )
   : scalars_( scalars )
   , forms_( forms )
   {
      WALBERLA_CHECK_EQUAL( scalars.size(), forms.size(), "Must pass same number of forms and scalars!" )
   }

   virtual ~N1E1LinearCombinationForm() = default;

   // --------------------------------------------------
   //  3D versions for tetrahedra (using new interface)
   // --------------------------------------------------

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix6r& elMat ) const override
   {
      elMat.setAll( 0 );
      for ( uint_t i = 0; i < forms_.size(); i++ )
      {
         Matrix6r tmpOut;
         forms_[i]->integrateAll( coords, tmpOut );
         elMat += scalars_[i] * tmpOut;
      }
   }

   bool assemble2D() const override
   {
      bool assemble = true;
      for ( const auto& form : forms_ )
      {
         assemble &= form->assemble2D();
      }
      return assemble;
   }

   bool assemble3D() const override
   {
      bool assemble = true;
      for ( const auto& form : forms_ )
      {
         assemble &= form->assemble3D();
      }
      return assemble;
   }

   bool assembly2DDefined() const override
   {
      bool assemble = true;
      for ( const auto& form : forms_ )
      {
         assemble &= form->assembly2DDefined();
      }
      return assemble;
   }

   bool assembly3DDefined() const override
   {
      bool assemble = true;
      for ( const auto& form : forms_ )
      {
         assemble &= form->assembly3DDefined();
      }
      return assemble;
   }

   virtual void setGeometryMap( const std::shared_ptr< GeometryMap >& geometryMap )
   {
      for ( auto& form : forms_ )
      {
         form->setGeometryMap( geometryMap );
      }
   }

 private:
   std::vector< real_t >    scalars_;
   std::vector< N1E1Form* > forms_;
};

} // namespace n1e1
} // namespace hyteg
