/*
 * Copyright (c) 2024 Michael Zikeli.
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

// Restrictions & Prolongations
#include "hyteg/gridtransferoperators/P1toP1LinearProlongation.hpp"
#include "hyteg/gridtransferoperators/P1toP1LinearRestriction.hpp"

namespace hyteg::mixedPrecisionMT {

enum class ProlongationRegister
{
   P1Linear
};
enum class RestrictionRegister
{
   P1Linear
};

template < ProlongationRegister Prolongation, RestrictionRegister Restriction, typename ValueType >
struct GridTransformationTrait;

template < typename ValueType >
struct GridTransformationTrait< ProlongationRegister::P1Linear, RestrictionRegister::P1Linear, ValueType >
{
   static std::shared_ptr< P1toP1LinearProlongation< ValueType > > prolongation()
   {
      return std::make_shared< P1toP1LinearProlongation< ValueType > >();
   }

   static std::shared_ptr< P1toP1LinearRestriction< ValueType > > restriction()
   {
      return std::make_shared< P1toP1LinearRestriction< ValueType > >();
   }
};

} // namespace hyteg::mixedPrecisionMT
