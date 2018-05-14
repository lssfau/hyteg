
#include "BoundaryConditions.hpp"
#include "core/debug/CheckFunctions.h"

namespace hhg {

BoundaryCondition BoundaryCondition::create012BC()
{
  BoundaryCondition bc;
  bc.createDirichletBC( "Dirichlet", 1 );
  bc.createNeumannBC( "Neumann", 2 );
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

DoFType BoundaryCondition::getBoundaryType( const uint_t & meshBoundaryFlag ) const
{
  if ( meshFlagToID_.count( meshBoundaryFlag ) == 0 )
  {
    return defaultBC_;
  }

  BoundaryUID boundaryUID = meshFlagToID_.at( meshBoundaryFlag );
  return boundaryUIDToType_.at( boundaryUID );
}

}