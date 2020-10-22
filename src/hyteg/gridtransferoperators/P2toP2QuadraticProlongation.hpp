/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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

#include "hyteg/gridtransferoperators/ProlongationOperator.hpp"

namespace hyteg {

template < typename VType >
class P2Function;

class P2toP2QuadraticProlongation : public ProlongationOperator< P2Function< walberla::real_t > >
{
 public:
   void prolongate( const P2Function< walberla::real_t >& function,
                    const walberla::uint_t&               sourceLevel,
                    const DoFType&                        flag ) const override;

   void prolongateAndAdd( const P2Function< walberla::real_t >& function,
                          const walberla::uint_t&               sourceLevel,
                          const DoFType&                        flag ) const override;

 private:
   void prolongateAdditively( const P2Function< walberla::real_t >& function,
                              const walberla::uint_t&               sourceLevel,
                              const DoFType&                        flag,
                              const UpdateType&                     updateType ) const;

   void prolongateAdditively3D( const P2Function< walberla::real_t >& function,
                                const walberla::uint_t&               sourceLevel,
                                const DoFType&                        flag,
                                const UpdateType&                     updateType ) const;

   void prolongateStandard( const P2Function< walberla::real_t >& function,
                            const walberla::uint_t&               sourceLevel,
                            const DoFType&                        flag ) const;
};

} // namespace hyteg
