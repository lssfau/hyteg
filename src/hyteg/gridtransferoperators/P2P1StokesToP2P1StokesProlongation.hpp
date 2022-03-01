/*
 * Copyright (c) 2017-2020 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/composites/P2P1TaylorHoodFunction.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P2toP2QuadraticProlongation.hpp"
#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"

namespace hyteg {

class P2P1StokesToP2P1StokesProlongation : public ProlongationOperator< P2P1TaylorHoodFunction< real_t > >
{
 public:
   typedef P2toP2QuadraticProlongation VelocityProlongation_T;
   typedef P1toP1LinearProlongation    PressureProlongation_T;

   void prolongate( const P2P1TaylorHoodFunction< real_t >& function,
                    const uint_t&                           sourceLevel,
                    const DoFType&                          flag ) const override
   {
      for ( uint_t k = 0; k < function.uvw().getDimension(); k++ )
      {
         quadraticProlongationOperator_.prolongate( function.uvw()[k], sourceLevel, flag );
      }
      linearProlongationOperator_.prolongate( function.p(), sourceLevel, flag );
   }

   void prolongateAndAdd( const P2P1TaylorHoodFunction< real_t >& function,
                          const uint_t&                           sourceLevel,
                          const DoFType&                          flag ) const override
   {
      for ( uint_t k = 0; k < function.uvw().getDimension(); k++ )
      {
         quadraticProlongationOperator_.prolongateAndAdd( function.uvw()[k], sourceLevel, flag );
      }
      linearProlongationOperator_.prolongateAndAdd( function.p(), sourceLevel, flag );
   }

 private:
   P2toP2QuadraticProlongation quadraticProlongationOperator_;
   P1toP1LinearProlongation    linearProlongationOperator_;
};
} // namespace hyteg
