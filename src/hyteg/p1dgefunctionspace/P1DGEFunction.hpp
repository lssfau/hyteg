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

#include "hyteg/dgfunctionspace/DGBasisLinearLagrange_Example.hpp"
#include "hyteg/dgfunctionspace/DGFunction.hpp"
#include "hyteg/p1functionspace/P1VectorFunction.hpp"
#include "hyteg/primitivestorage/PrimitiveStorage.hpp"

namespace hyteg {

template < typename ValueType >
class P1DGEFunction final : public Function< P1DGEFunction< ValueType > >
{
 public:
   P1DGEFunction( const std::string&                         name,
                  const std::shared_ptr< PrimitiveStorage >& storage,
                  uint_t                                     minLevel,
                  uint_t                                     maxLevel,
                  BoundaryCondition                          boundaryCondition = BoundaryCondition::create0123BC() );

   std::shared_ptr< P1VectorFunction< ValueType > > getConformingPart() const { return u_conforming_; }

   std::shared_ptr< dg::DGFunction< ValueType > > getDiscontinuousPart() const { return u_discontinuous_; }

 protected:
   std::shared_ptr< dg::DGBasisLinearLagrange_Example > basis_;

   std::shared_ptr< P1VectorFunction< ValueType > > u_conforming_;
   std::shared_ptr< dg::DGFunction< ValueType > >   u_discontinuous_;
};

template < typename ValueType >
P1DGEFunction< ValueType >::P1DGEFunction( const std::string&                         name,
                                           const std::shared_ptr< PrimitiveStorage >& storage,
                                           uint_t                                     minLevel,
                                           uint_t                                     maxLevel,
                                           BoundaryCondition                          boundaryCondition )
: Function< P1DGEFunction< ValueType > >( name, storage, minLevel, maxLevel )
, basis_{ std::make_shared< dg::DGBasisLinearLagrange_Example >() }
, u_conforming_{ std::make_shared< P1VectorFunction< ValueType > >( name, storage, minLevel, maxLevel, boundaryCondition ) }
, u_discontinuous_{
      std::make_shared< dg::DGFunction< ValueType > >( name, storage, minLevel, maxLevel, basis_, 0, boundaryCondition ) }
{}

} // namespace hyteg
