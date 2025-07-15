/*
 * Copyright (c) 2022 Marcus Mohr.
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

#include "hyteg/fenics/fenics.hpp"
#include "hyteg/fenics/ufc_traits.hpp"
#include "hyteg/forms/P1Form.hpp"

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

// N1curl
#include "hyteg/forms/form_fenics_generated/n1Curl_tet_mass.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

namespace hyteg {

template < class UFCOperator2D, class UFCOperator3D = fenics::UndefinedAssembly >
class N1curlFenicsForm : public Form
{
 public:

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix6r& elMat ) const
   {
      double fenicsCoords[12];
      fenicsCoords[0] = coords[0][0];
      fenicsCoords[1] = coords[0][1];
      fenicsCoords[2] = coords[0][2];

      fenicsCoords[3] = coords[1][0];
      fenicsCoords[4] = coords[1][1];
      fenicsCoords[5] = coords[1][2];

      fenicsCoords[6] = coords[2][0];
      fenicsCoords[7] = coords[2][1];
      fenicsCoords[8] = coords[2][2];

      fenicsCoords[9]  = coords[3][0];
      fenicsCoords[10] = coords[3][1];
      fenicsCoords[11] = coords[3][2];

      UFCOperator3D gen;
      hyteg::Matrix< double, 6, 6 > tmp = elMat.cast< double >();
      gen.tabulate_tensor( tmp.data(), nullptr, fenicsCoords, 0 );
      elMat = tmp.cast< real_t >();
   }

   void setGeometryMap( const std::shared_ptr< GeometryMap >& ) const override {}
};

} // namespace hyteg
