/*
 * Copyright (c) 2017-2020 Nils Kohl.
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

#include "hyteg/petsc/PETScWrapper.hpp"
#include "hyteg/sparseassembly/VectorProxy.hpp"

namespace hyteg {

using walberla::uint_t;

class PETScVectorProxy : public VectorProxy
{
 public:
   PETScVectorProxy( Vec vec )
   : vec_( vec )
   {}

   /// \brief Sets the passed value in the vector.
   virtual void setValue( uint_t idx, real_t value )
   {
      VecSetValue( vec_, static_cast< PetscInt >( idx ), value, INSERT_VALUES );
   }

   /// \brief Returns the passed value of the vector.
   virtual real_t getValue( uint_t idx ) const
   {
      PetscInt idxPetsc = static_cast< PetscInt >( idx );
      real_t   value;
      VecGetValues( vec_, 1, &idxPetsc, &value );
      return value;
   }

 private:
   Vec vec_;
};

} // namespace hyteg