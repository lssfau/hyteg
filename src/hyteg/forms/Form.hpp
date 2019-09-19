/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Marcus Mohr.
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

#include "hyteg/geometry/GeometryMap.hpp"

namespace hyteg {

class Form
{
 public:
   virtual ~Form() {}

   virtual bool assemble2D() const = 0;

   virtual bool assemble3D() const = 0;

   virtual bool assembly2DDefined() const = 0;

   virtual bool assembly3DDefined() const = 0;

   std::shared_ptr< GeometryMap > geometryMap;
};

} // namespace hyteg