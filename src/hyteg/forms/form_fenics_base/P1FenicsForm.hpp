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

#include "hyteg/fenics/fenics.hpp"
#include "hyteg/fenics/ufc_traits.hpp"
#include "hyteg/forms/P1Form.hpp"

#ifdef _MSC_VER
#pragma warning( push, 0 )
#endif

// P1
#include "hyteg/forms/form_fenics_generated/p1_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p1_div.h"
#include "hyteg/forms/form_fenics_generated/p1_div_K_grad.h"
#include "hyteg/forms/form_fenics_generated/p1_divt.h"
#include "hyteg/forms/form_fenics_generated/p1_mass.h"
#include "hyteg/forms/form_fenics_generated/p1_polar_laplacian.h"
#include "hyteg/forms/form_fenics_generated/p1_polar_mass.h"
#include "hyteg/forms/form_fenics_generated/p1_pspg.h"
#include "hyteg/forms/form_fenics_generated/p1_stokes_epsilon.h"
#include "hyteg/forms/form_fenics_generated/p1_stokes_full.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_div_tet.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_divt_tet.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_mass.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_pspg_tet.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_stokes_epsilon_tet.h"
#include "hyteg/forms/form_fenics_generated/p1_tet_stokes_full_tet.h"

// P2
#include "hyteg/forms/form_fenics_generated/p2_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p2_div.h"
#include "hyteg/forms/form_fenics_generated/p2_divt.h"
#include "hyteg/forms/form_fenics_generated/p2_mass.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_diffusion.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_div_tet.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_divt_tet.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_mass.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_pspg_tet.h"
#include "hyteg/forms/form_fenics_generated/p2_tet_stokes_epsilon_tet.h"

// P1 to P2
#include "hyteg/forms/form_fenics_generated/p1_to_p2_divt.h"
#include "hyteg/forms/form_fenics_generated/p1_to_p2_tet_divt_tet.h"

// P2 to P1
#include "hyteg/forms/form_fenics_generated/p2_to_p1_div.h"
#include "hyteg/forms/form_fenics_generated/p2_to_p1_tet_div_tet.h"

#ifdef _MSC_VER
#pragma warning( pop )
#endif

namespace hyteg {

template < class UFCOperator2D, class UFCOperator3D = fenics::UndefinedAssembly >
class P1FenicsForm : public P1Form
{
 public:
   void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const override
   {
      hyteg::Matrix< double, 3, 3 > localStiffnessMatrix{ hyteg::Matrix< double, 3, 3 >::Zero() };
      double                        fenicsCoords[6];
      fenicsCoords[0] = coords[0][0];
      fenicsCoords[1] = coords[0][1];
      fenicsCoords[2] = coords[1][0];
      fenicsCoords[3] = coords[1][1];
      fenicsCoords[4] = coords[2][0];
      fenicsCoords[5] = coords[2][1];
      UFCOperator2D gen;
      gen.tabulate_tensor( localStiffnessMatrix.data(), nullptr, fenicsCoords, 0 );
      out[0] = localStiffnessMatrix( 0, 0 );
      out[1] = localStiffnessMatrix( 0, 1 );
      out[2] = localStiffnessMatrix( 0, 2 );
   }

   void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const override
   {
      //typename fenics::UFCTrait< UFCOperator3D >::LocalStiffnessMatrix_T localStiffnessMatrix;

      // Flattening the offset array to be able to pass it to the fenics routines.
      double geometricOffsetsArray[12];
      for ( uint_t cellVertex = 0; cellVertex < 4; cellVertex++ )
      {
         for ( int coordinate = 0; coordinate < 3; coordinate++ )
         {
            geometricOffsetsArray[cellVertex * 3 + walberla::uint_c( coordinate )] = coords[cellVertex][coordinate];
         }
      }

      UFCOperator3D gen;
      Eigen::Matrix< double,
                     fenics::UFCTrait< UFCOperator3D >::LocalStiffnessMatrix_T::RowsAtCompileTime,
                     fenics::UFCTrait< UFCOperator3D >::LocalStiffnessMatrix_T::ColsAtCompileTime,
                     Eigen::RowMajor >
          localStiffnessMatrix;
      gen.tabulate_tensor( localStiffnessMatrix.data(), NULL, geometricOffsetsArray, 0 );

      out[0] = localStiffnessMatrix( 0, 0 );
      out[1] = localStiffnessMatrix( 0, 1 );
      out[2] = localStiffnessMatrix( 0, 2 );
      out[3] = localStiffnessMatrix( 0, 3 );
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

   void integrateRow0( const std::array< Point3D, 4 >& coords, Matrixr< 1, 4 >& elMat ) const override
   {
      Point4D row;
      integrate(coords, row);
      for (int i = 0; i < 4; ++i)
      {
         elMat(0,i) = row[i];
      }
   }

   void integrateAll( const std::array< Point3D, 3 >& coords, Matrix3r& elMat ) const override
   {
      double fenicsCoords[6];
      fenicsCoords[0] = coords[0][0];
      fenicsCoords[1] = coords[0][1];
      fenicsCoords[2] = coords[1][0];
      fenicsCoords[3] = coords[1][1];
      fenicsCoords[4] = coords[2][0];
      fenicsCoords[5] = coords[2][1];
      UFCOperator2D                 gen;
      hyteg::Matrix< double, 3, 3 > tmp = elMat.cast< double >();
      gen.tabulate_tensor( tmp.data(), nullptr, fenicsCoords, 0 );
      elMat = tmp.cast< real_t >();
   }

   void integrateAll( const std::array< Point3D, 4 >& coords, Matrix4r& elMat ) const override
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

      UFCOperator3D                 gen;
      hyteg::Matrix< double, 4, 4 > tmp = elMat.cast< double >();
      gen.tabulate_tensor( tmp.data(), nullptr, fenicsCoords, 0 );
      elMat = tmp.cast< real_t >();
   }

   inline void setGeometryMap( const std::shared_ptr< GeometryMap >& map ) const { WALBERLA_UNUSED( map ); }
};

} // namespace hyteg
