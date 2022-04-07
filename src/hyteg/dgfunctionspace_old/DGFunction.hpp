/*
 * Copyright (c) 2017-2019 Boerge Struempfel, Daniel Drzisga, Dominik Thoennes, Nils Kohl.
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

#include "hyteg/facedofspace_old/FaceDoFFunction.hpp"

namespace hyteg {

template < typename ValueType >
class DGFunction_old : public FaceDoFFunction_old< ValueType >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = DGFunction_old< VType >;

   DGFunction_old( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel )
   : FaceDoFFunction_old< ValueType >( name, storage, minLevel, maxLevel )
   {}

   void projectP1( P1Function< ValueType >& src, uint_t level, DoFType flag, UpdateType updateType = Replace );
};

} // namespace hyteg
