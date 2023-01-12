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

#include "hyteg/forms/P2Form.hpp"

namespace hyteg {

/// \brief P2 form that performs a linear combination of other P2 forms.
///
/// The linear combinations are formed according to the passed scalars and forms:
///    result = scalar[0] * form[0]->integrate() + scalar[1] * form[1]->integrate() + ...
/// The number of passed scalars and forms must be identical.
///
class P2LinearCombinationForm : public P2Form
{
 public:
   P2LinearCombinationForm() = default;

   P2LinearCombinationForm( const std::vector< real_t >& scalars, const std::vector< P2Form* >& forms )
   : scalars_( scalars )
   , forms_( forms )
   {
      WALBERLA_CHECK_EQUAL( scalars.size(), forms.size(), "Must pass same number of forms and scalars!" )
   }

   virtual ~P2LinearCombinationForm() = default;

   // ---------------------------
   //  2D versions for triangles
   // ---------------------------
   void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      out.setAll( 0 );
      for ( uint_t i = 0; i < forms_.size(); i++ )
      {
         Point3D tmpOut;
         forms_[i]->integrate( coords, tmpOut );
         out += scalars_[i] * tmpOut;
      }
   }

   void integrateEdgeToVertex( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      out.setAll( 0 );
      for ( uint_t i = 0; i < forms_.size(); i++ )
      {
         Point3D tmpOut;
         forms_[i]->integrateEdgeToVertex( coords, tmpOut );
         out += scalars_[i] * tmpOut;
      }
   }

   void integrateVertexToEdge( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      out.setAll( 0 );
      for ( uint_t i = 0; i < forms_.size(); i++ )
      {
         Point3D tmpOut;
         forms_[i]->integrateVertexToEdge( coords, tmpOut );
         out += scalars_[i] * tmpOut;
      }
   }

   void integrateEdgeToEdge( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      out.setAll( 0 );
      for ( uint_t i = 0; i < forms_.size(); i++ )
      {
         Point3D tmpOut;
         forms_[i]->integrateEdgeToEdge( coords, tmpOut );
         out += scalars_[i] * tmpOut;
      }
   }

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix6r& elMat ) const
   {
      elMat.setAll( 0 );
      for ( uint_t i = 0; i < forms_.size(); i++ )
      {
         Matrix6r tmpOut;
         forms_[i]->integrateAll( coords, tmpOut );
         elMat += scalars_[i] * tmpOut;
      }
   }

   // ----------------------------
   //  3D versions for tetrahedra
   // ----------------------------
   void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      out.setAll( 0 );
      for ( uint_t i = 0; i < forms_.size(); i++ )
      {
         Point4D tmpOut;
         forms_[i]->integrate( coords, tmpOut );
         out += scalars_[i] * tmpOut;
      }
   }

   void integrateEdgeToVertex( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2LinearCombinationFenicsForm not implemented for 3D!" );
   }

   void integrateVertexToEdge( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2LinearCombinationFenicsForm not implemented for 3D!" );
   }

   void integrateEdgeToEdge( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2LinearCombinationFenicsForm not implemented for 3D!" );
   }

   // --------------------------------------------------
   //  3D versions for tetrahedra (using new interface)
   // --------------------------------------------------

   real_t integrate( const std::array< Point3D, 4 >&     coords,
                     const P2Form::dofPosByVertexPair3D& cntrPos,
                     const P2Form::dofPosByVertexPair3D& leafPos ) const
   {
      real_t out = 0;
      for ( uint_t i = 0; i < forms_.size(); i++ )
      {
         out += scalars_[i] * forms_[i]->integrate( coords, cntrPos, leafPos );
      }
      return out;
   }

   std::vector< real_t > integrate( const std::array< Point3D, 4 >&                    coords,
                                    const P2Form::dofPosByVertexPair3D&                cntrPos,
                                    const std::vector< P2Form::dofPosByVertexPair3D >& leafPos ) const
   {
      std::vector< real_t > out( leafPos.size(), 0 );
      for ( uint_t i = 0; i < forms_.size(); i++ )
      {
         auto tmpOut = forms_[i]->integrate( coords, cntrPos, leafPos );
         for ( uint_t j = 0; j < leafPos.size(); j++ )
         {
            out[j] += scalars_[i] * tmpOut[j];
         }
      }
      return out;
   }

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix10r& elMat ) const
   {
      elMat.setAll( 0 );
      for ( uint_t i = 0; i < forms_.size(); i++ )
      {
         Matrix10r tmpOut;
         forms_[i]->integrateAll( coords, tmpOut );
         elMat += scalars_[i] * tmpOut;
      }
   }

   virtual void setGeometryMap( const std::shared_ptr< GeometryMap >& geometryMap )
   {
      for ( auto& form : forms_ )
      {
         form->setGeometryMap( geometryMap );
      }
   }

 private:
   std::vector< real_t >  scalars_;
   std::vector< P2Form* > forms_;
};

} // namespace hyteg
