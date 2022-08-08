/*
* Copyright (c) 2017-2022 Nils Kohl.
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

#include "hyteg/types/pointnd.hpp"

#include "Eigen/Eigen"

namespace hyteg {

class EGBasis {
 public:

   void integrateBasisFunction( uint_t                                                degree,
                                const std::array< Eigen::Matrix< real_t, 2, 1 >, 3 >& coords,
                                const std::function< real_t( const Point3D& ) >&      f0,
                                const std::function< real_t( const Point3D& ) >&      f1,
                                std::vector< real_t >&                                values )
   {
      WALBERLA_CHECK(degree == 0);
      WALBERLA_CHECK(values.size() == 1);

      const auto p_affine_0_0 = coords[0]( 0 );
      const auto p_affine_0_1 = coords[0]( 1 );

      const auto p_affine_1_0 = coords[1]( 0 );
      const auto p_affine_1_1 = coords[1]( 1 );

      const auto p_affine_2_0 = coords[2]( 0 );
      const auto p_affine_2_1 = coords[2]( 1 );

      callback_Scalar_Variable_Coefficient_2D_f0 = f0;
      callback_Scalar_Variable_Coefficient_2D_f1 = f1;

      real_t Scalar_Variable_Coefficient_2D_f0_out0_id0 = 0;
      real_t Scalar_Variable_Coefficient_2D_f1_out0_id1 = 0;
      real_t Scalar_Variable_Coefficient_2D_f0_out0_id2 = 0;
      real_t Scalar_Variable_Coefficient_2D_f1_out0_id3 = 0;
      real_t Scalar_Variable_Coefficient_2D_f0_out0_id4 = 0;
      real_t Scalar_Variable_Coefficient_2D_f1_out0_id5 = 0;
      real_t Scalar_Variable_Coefficient_2D_f0_out0_id6 = 0;
      real_t Scalar_Variable_Coefficient_2D_f1_out0_id7 = 0;
      real_t Scalar_Variable_Coefficient_2D_f0_out0_id8 = 0;
      real_t Scalar_Variable_Coefficient_2D_f1_out0_id9 = 0;
      real_t Scalar_Variable_Coefficient_2D_f0_out0_id10 = 0;
      real_t Scalar_Variable_Coefficient_2D_f1_out0_id11 = 0;

      Scalar_Variable_Coefficient_2D_f0( 0.091576213509770743*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.81684757298045851*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.81684757298045851*p_affine_2_1, &Scalar_Variable_Coefficient_2D_f0_out0_id0 );
      Scalar_Variable_Coefficient_2D_f1( 0.091576213509770743*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.81684757298045851*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.81684757298045851*p_affine_2_1, &Scalar_Variable_Coefficient_2D_f1_out0_id1 );
      Scalar_Variable_Coefficient_2D_f0( 0.44594849091596489*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.10810301816807022*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.10810301816807022*p_affine_2_1, &Scalar_Variable_Coefficient_2D_f0_out0_id2 );
      Scalar_Variable_Coefficient_2D_f1( 0.44594849091596489*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.10810301816807022*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.10810301816807022*p_affine_2_1, &Scalar_Variable_Coefficient_2D_f1_out0_id3 );
      Scalar_Variable_Coefficient_2D_f0( 0.091576213509770743*p_affine_0_0 + 0.81684757298045851*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.81684757298045851*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_f0_out0_id4 );
      Scalar_Variable_Coefficient_2D_f1( 0.091576213509770743*p_affine_0_0 + 0.81684757298045851*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.091576213509770743*p_affine_0_1 + 0.81684757298045851*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_f1_out0_id5 );
      Scalar_Variable_Coefficient_2D_f0( 0.44594849091596489*p_affine_0_0 + 0.10810301816807022*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.10810301816807022*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_f0_out0_id6 );
      Scalar_Variable_Coefficient_2D_f1( 0.44594849091596489*p_affine_0_0 + 0.10810301816807022*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.44594849091596489*p_affine_0_1 + 0.10810301816807022*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_f1_out0_id7 );
      Scalar_Variable_Coefficient_2D_f0( 0.81684757298045851*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.81684757298045851*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_f0_out0_id8 );
      Scalar_Variable_Coefficient_2D_f1( 0.81684757298045851*p_affine_0_0 + 0.091576213509770743*p_affine_1_0 + 0.091576213509770743*p_affine_2_0, 0.81684757298045851*p_affine_0_1 + 0.091576213509770743*p_affine_1_1 + 0.091576213509770743*p_affine_2_1, &Scalar_Variable_Coefficient_2D_f1_out0_id9 );
      Scalar_Variable_Coefficient_2D_f0( 0.10810301816807022*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.10810301816807022*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_f0_out0_id10 );
      Scalar_Variable_Coefficient_2D_f1( 0.10810301816807022*p_affine_0_0 + 0.44594849091596489*p_affine_1_0 + 0.44594849091596489*p_affine_2_0, 0.10810301816807022*p_affine_0_1 + 0.44594849091596489*p_affine_1_1 + 0.44594849091596489*p_affine_2_1, &Scalar_Variable_Coefficient_2D_f1_out0_id11 );

      real_t tmp_0 = std::abs(p_affine_0_0*p_affine_1_1 - p_affine_0_0*p_affine_2_1 - p_affine_0_1*p_affine_1_0 + p_affine_0_1*p_affine_2_0 + p_affine_1_0*p_affine_2_1 - p_affine_1_1*p_affine_2_0);
      real_t a_0_0 = 0.054975871827660928*tmp_0*(-0.24175711982356257*Scalar_Variable_Coefficient_2D_f0_out0_id0 + 0.4835142396471252*Scalar_Variable_Coefficient_2D_f1_out0_id1) + 0.11169079483900572*tmp_0*(0.11261515758263158*Scalar_Variable_Coefficient_2D_f0_out0_id10 + 0.11261515758263158*Scalar_Variable_Coefficient_2D_f1_out0_id11) + 0.11169079483900572*tmp_0*(0.11261515758263158*Scalar_Variable_Coefficient_2D_f0_out0_id2 - 0.2252303151652631*Scalar_Variable_Coefficient_2D_f1_out0_id3) + 0.054975871827660928*tmp_0*(0.4835142396471252*Scalar_Variable_Coefficient_2D_f0_out0_id4 - 0.24175711982356257*Scalar_Variable_Coefficient_2D_f1_out0_id5) + 0.11169079483900572*tmp_0*(-0.2252303151652631*Scalar_Variable_Coefficient_2D_f0_out0_id6 + 0.11261515758263158*Scalar_Variable_Coefficient_2D_f1_out0_id7) + 0.054975871827660928*tmp_0*(-0.24175711982356257*Scalar_Variable_Coefficient_2D_f0_out0_id8 - 0.24175711982356257*Scalar_Variable_Coefficient_2D_f1_out0_id9);

      values[0] = a_0_0;
   }

 private:
   void Scalar_Variable_Coefficient_2D_f0( real_t in_0, real_t in_1, real_t* out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_f0( Point3D( { in_0, in_1, 0 } ) );
   }

   void Scalar_Variable_Coefficient_2D_f1( real_t in_0, real_t in_1, real_t* out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_f1( Point3D( { in_0, in_1, 0 } ) );
   }

   void Scalar_Variable_Coefficient_2D_f2( real_t in_0, real_t in_1, real_t* out_0 ) const
   {
      *out_0 = callback_Scalar_Variable_Coefficient_2D_f2( Point3D( { in_0, in_1, 0 } ) );
   }

   std::function< real_t( const Point3D& ) > callback_Scalar_Variable_Coefficient_2D_f0;
   std::function< real_t( const Point3D& ) > callback_Scalar_Variable_Coefficient_2D_f1;
   std::function< real_t( const Point3D& ) > callback_Scalar_Variable_Coefficient_2D_f2;
};

}