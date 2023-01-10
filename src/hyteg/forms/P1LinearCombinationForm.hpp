/*
 * Copyright (c) 2017-2020 Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/forms/P1Form.hpp"

namespace hyteg {

/// \brief P1 form that performs a linear combination of other P1 forms.
///
/// The linear combinations are formed according to the passed scalars and forms:
///    result = scalar[0] * form[0]->integrate() + scalar[1] * form[1]->integrate() + ...
/// The number of passed scalars and forms must be identical.
///
class P1LinearCombinationForm : public P1Form
{
 public:
   P1LinearCombinationForm() = default;

   P1LinearCombinationForm( const std::vector< real_t >& scalars, const std::vector< P1Form* >& forms )
       : scalars_unowned_forms_( scalars )
       , unowned_forms_( forms )
   {
      WALBERLA_CHECK_EQUAL( scalars.size(), forms.size(), "Must pass same number of forms and scalars!" )
   }

   virtual ~P1LinearCombinationForm() = default;

   template < typename NewFormType, typename... ConstructorArguments >
   void addOwnedForm( const real_t& scalar, ConstructorArguments... args )
   {
      scalars_owned_forms_.push_back( scalar );
      owned_forms_.push_back( std::make_shared< NewFormType >( args... ) );
   }

   void addUnownedForm( const real_t& scalar, P1Form* form )
   {
      scalars_unowned_forms_.push_back( scalar );
      unowned_forms_.push_back( form );
   }

   // ---------------------------
   //  2D versions for triangles
   // ---------------------------
   void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      out.setAll( 0 );
      for ( uint_t i = 0; i < unowned_forms_.size(); i++ )
      {
         Point3D tmpOut;
         unowned_forms_[i]->integrate( coords, tmpOut );
         out += scalars_unowned_forms_[i] * tmpOut;
      }
      for ( uint_t i = 0; i < owned_forms_.size(); i++ )
      {
         Point3D tmpOut;
         owned_forms_[i]->integrate( coords, tmpOut );
         out += scalars_owned_forms_[i] * tmpOut;
      }
   }

   void integrateRow0( const std::array< Point3D, 3 >& coords, Matrixr< 1, 3 >& elMat ) const override
   {
      Point3D row;
      integrate(coords, row);
      for (int i = 0; i < 3; ++i)
      {
         elMat(0,i) = row[i];
      }
   }

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix3r& elMat ) const
   {
      elMat.setAll( 0 );
      for ( uint_t i = 0; i < unowned_forms_.size(); i++ )
      {
         Matrix3r tmpOut;
         unowned_forms_[i]->integrateAll( coords, tmpOut );
         elMat += scalars_unowned_forms_[i] * tmpOut;
      }
      for ( uint_t i = 0; i < owned_forms_.size(); i++ )
      {
         Matrix3r tmpOut;
         owned_forms_[i]->integrateAll( coords, tmpOut );
         elMat += scalars_owned_forms_[i] * tmpOut;
      }
   }

   // ----------------------------
   //  3D versions for tetrahedra
   // ----------------------------
   void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      out.setAll( 0 );
      for ( uint_t i = 0; i < unowned_forms_.size(); i++ )
      {
         Point4D tmpOut;
         unowned_forms_[i]->integrate( coords, tmpOut );
         out += scalars_unowned_forms_[i] * tmpOut;
      }
      for ( uint_t i = 0; i < owned_forms_.size(); i++ )
      {
         Point4D tmpOut;
         owned_forms_[i]->integrate( coords, tmpOut );
         out += scalars_owned_forms_[i] * tmpOut;
      }
   }

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrixr< 1, 4 >& elMat ) const override
   {
      Point4D row;
      integrate(coords, row);
      for (int i = 0; i < 4; ++i)
      {
         elMat(0,i) = row[i];
      }
   }

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix4r& elMat ) const
   {
      elMat.setAll( 0 );
      for ( uint_t i = 0; i < unowned_forms_.size(); i++ )
      {
         Matrix4r tmpOut;
         unowned_forms_[i]->integrateAll( coords, tmpOut );
         elMat += scalars_unowned_forms_[i] * tmpOut;
      }
      for ( uint_t i = 0; i < owned_forms_.size(); i++ )
      {
         Matrix4r tmpOut;
         owned_forms_[i]->integrateAll( coords, tmpOut );
         elMat += scalars_owned_forms_[i] * tmpOut;
      }
   }

   virtual void setGeometryMap( const std::shared_ptr< GeometryMap >& geometryMap )
   {
      for ( auto& form : unowned_forms_ )
      {
         form->setGeometryMap( geometryMap );
      }
      for ( auto& form : owned_forms_ )
      {
         form->setGeometryMap( geometryMap );
      }
   }

 private:
   std::vector< real_t >  scalars_unowned_forms_;
   std::vector< P1Form* > unowned_forms_;

   std::vector< real_t >                    scalars_owned_forms_;
   std::vector< std::shared_ptr< P1Form > > owned_forms_;
};

} // namespace hyteg
