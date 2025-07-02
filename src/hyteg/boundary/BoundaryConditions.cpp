/*
 * Copyright (c) 2017-2025 Dominik Thoennes, Nils Kohl, Andreas Burkhart.
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

#include "BoundaryConditions.hpp"
#include "core/debug/CheckFunctions.h"

namespace hyteg {

BoundaryCondition BoundaryCondition::create0123BC()
{
  BoundaryCondition bc;
  bc.createDirichletBC( "Dirichlet", 1 );
  bc.createNeumannBC( "Neumann", 2 );
  bc.createFreeslipBC( "Freeslip", 3 );
  return bc;
}

BoundaryCondition BoundaryCondition::createAllInnerBC()
{
  BoundaryCondition bc( hyteg::Inner );
  return bc;
}

BoundaryUID BoundaryCondition::createDirichletBC( const std::string & name, const uint_t & meshBoundaryFlag )
{
  std::vector< uint_t > meshBoundaryFlags = { meshBoundaryFlag };
  return createDirichletBC( name, meshBoundaryFlags );
}

BoundaryUID BoundaryCondition::createDirichletBC( const std::string & name, const std::vector< uint_t > & meshBoundaryFlags )
{
  BoundaryUID boundaryUID( name );
  WALBERLA_CHECK_EQUAL( boundaryUIDToType_.count( boundaryUID ), 0, "Boundary condition with name '" << name << "' already exists!" )
  boundaryUIDToType_[ boundaryUID ] = DoFType::DirichletBoundary;

  for ( const auto & meshBoundaryFlag : meshBoundaryFlags )
  {
    meshFlagToID_[ meshBoundaryFlag ] = boundaryUID;
  }

  return boundaryUID;
}

BoundaryUID BoundaryCondition::createNeumannBC( const std::string & name, const uint_t & meshBoundaryFlag )
{
  std::vector< uint_t > meshBoundaryFlags = { meshBoundaryFlag };
  return createNeumannBC( name, meshBoundaryFlags );
}

BoundaryUID BoundaryCondition::createNeumannBC( const std::string & name, const std::vector< uint_t > & meshBoundaryFlags )
{
  BoundaryUID boundaryUID( name );
  WALBERLA_CHECK_EQUAL( boundaryUIDToType_.count( boundaryUID ), 0, "Boundary condition with name '" << name << "' already exists!" )
  boundaryUIDToType_[ boundaryUID ] = DoFType::NeumannBoundary;

  for ( const auto & meshBoundaryFlag : meshBoundaryFlags )
  {
    meshFlagToID_[ meshBoundaryFlag ] = boundaryUID;
  }

  return boundaryUID;
}

BoundaryUID BoundaryCondition::createFreeslipBC( const std::string & name, const uint_t & meshBoundaryFlag )
{
   std::vector< uint_t > meshBoundaryFlags = { meshBoundaryFlag };
   return createFreeslipBC( name, meshBoundaryFlags );
}

BoundaryUID BoundaryCondition::createFreeslipBC( const std::string & name, const std::vector< uint_t > & meshBoundaryFlags )
{
   BoundaryUID boundaryUID( name );
   WALBERLA_CHECK_EQUAL( boundaryUIDToType_.count( boundaryUID ), 0, "Boundary condition with name '" << name << "' already exists!" )
   boundaryUIDToType_[ boundaryUID ] = DoFType::FreeslipBoundary;

   for ( const auto & meshBoundaryFlag : meshBoundaryFlags )
   {
      meshFlagToID_[ meshBoundaryFlag ] = boundaryUID;
   }

   return boundaryUID;
}

BoundaryUID BoundaryCondition::createInnerBC( const std::string & name, const uint_t &meshBoundaryFlag )
{
   std::vector< uint_t > meshBoundaryFlags = { meshBoundaryFlag };
   return createInnerBC( name, meshBoundaryFlags );
}

BoundaryUID BoundaryCondition::createInnerBC( const std::string & name, const std::vector< uint_t > & meshBoundaryFlags )
{
   BoundaryUID boundaryUID( name );
   WALBERLA_CHECK_EQUAL( boundaryUIDToType_.count( boundaryUID ), 0, "Boundary condition with name '" << name << "' already exists!" )
   boundaryUIDToType_[ boundaryUID ] = DoFType::Inner;

   for ( const auto & meshBoundaryFlag : meshBoundaryFlags )
   {
      meshFlagToID_[ meshBoundaryFlag ] = boundaryUID;
   }

   return boundaryUID;
}

DoFType BoundaryCondition::getBoundaryType( const uint_t & meshBoundaryFlag ) const
{
  if ( meshFlagToID_.count( meshBoundaryFlag ) == 0 )
  {
    return defaultBC_;
  }

  BoundaryUID boundaryUID = meshFlagToID_.at( meshBoundaryFlag );
  return boundaryUIDToType_.at( boundaryUID );
}

bool BoundaryCondition::operator==( const hyteg::BoundaryCondition & other ) const
{
  return    defaultBC_         == other.defaultBC_
         && meshFlagToID_      == other.meshFlagToID_
         && boundaryUIDToType_ == other.boundaryUIDToType_;
}

}