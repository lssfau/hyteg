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

#include "hyteg/forms/N1E1Form.hpp"

namespace hyteg {

class N1E1FormHyTeG : public n1e1::N1E1Form
{
 public:
   virtual ~N1E1FormHyTeG() {}

   bool assemble2D() const override
   {
      WALBERLA_ABORT( "Don't call assemble2D on an N1E1FormHyteG child" );
      return false;
   };
   bool assemble3D() const override
   {
      WALBERLA_ABORT( "Don't call assemble3D on an N1E1FormHyteG child" );
      return false;
   };
   bool assembly2DDefined() const override
   {
      WALBERLA_ABORT( "Don't call assembly2DDefined on an N1E1FormHyteG child" );
      return false;
   };
   bool assembly3DDefined() const override
   {
      WALBERLA_ABORT( "Don't call assembly2DDefined on an N1E1FormHyteG child" );
      return false;
   };
};

} // namespace hyteg
