#pragma once

#include "core/DataTypes.h"
#include "core/uid/all.h"
#include "hyteg/types/flags.hpp"

namespace hyteg {

using walberla::uint_t;

class BoundaryUIDGenerator : public walberla::uid::IndexGenerator< BoundaryUIDGenerator, uint_t >{};
typedef walberla::UID< BoundaryUIDGenerator > BoundaryUID;

/// \brief A \ref BoundaryCondition maps mesh boundary flags to boundary conditions.
///
/// Mesh boundary flags are integers that are defined by a mesh file of generated \ref MeshInfo.
/// Each \ref BoundaryCondition instance represents a mapping of the mesh boundary flags to one set of boundary conditions (BCs).
/// This way, multiple sets of BCs can be defined in the same application and attached FE functions individually.
///
/// The class introduces the \ref BoundaryUID which is assigned to a set of mesh boundary flags that is set to a BC.
/// Through this UID, multiple sets of the same BC type (i.e. Dirichlet) can be distinguished.
class BoundaryCondition
{
public:

  /// Creates a \ref BoundaryCondition that assigns:
  ///
  /// mesh boundary flag | BC
  /// -------------------|-------
  ///                  1 | Dirichlet
  ///                  2 | Neumann
  ///              other | Domain
  ///
  static BoundaryCondition create012BC();
  static BoundaryCondition createAllInnerBC();

  /// Creates a \ref BoundaryCondition.
  /// \param defaultBC BC that is set to all mesh boundary flags that are not explicitly registered
  BoundaryCondition( const DoFType & defaultBC = DoFType::Inner ) :
    defaultBC_( defaultBC ), defaultBoundaryUID_( "DefaultBoundaryUID" )
  {}

  /// Return the BoundaryUID of the passed mesh flag.
  /// If no boundary UID was assigned, the default BoundaryUID is returned.
  BoundaryUID getBoundaryUIDFromMeshFlag( const uint_t & meshFlag ) const
  {
    if ( meshFlagToID_.count( meshFlag ) == 0 )
    {
      return defaultBoundaryUID_;
    }
    return meshFlagToID_.at( meshFlag );
  }

  /// Assigns one mesh boundary flag to the domain.
  BoundaryUID createDomainBC   ( const std::string & name, const uint_t &                meshBoundaryFlag  );
  /// Assigns multiple mesh boundary flags to the domain.
  BoundaryUID createDomainBC   ( const std::string & name, const std::vector< uint_t > & meshBoundaryFlags );

  /// Assigns one mesh boundary flag to a Dirichlet BC.
  BoundaryUID createDirichletBC( const std::string & name, const uint_t &                meshBoundaryFlag  );
  /// Assigns multiple mesh boundary flags to a Dirichlet BC.
  BoundaryUID createDirichletBC( const std::string & name, const std::vector< uint_t > & meshBoundaryFlags );

  /// Assigns one mesh boundary flag to a Neumann BC.
  BoundaryUID createNeumannBC  ( const std::string & name, const uint_t &                meshBoundaryFlag  );
  /// Assigns multiple mesh boundary flags to a Neumann BC.
  BoundaryUID createNeumannBC  ( const std::string & name, const std::vector< uint_t > & meshBoundaryFlags );

  /// Returns the boundary type that is assigned to the passed mesh boundary flag integer.
  DoFType getBoundaryType( const uint_t & meshBoundaryFlag ) const;

  bool operator== ( const BoundaryCondition & other ) const;

private:

  DoFType                          defaultBC_;
  BoundaryUID                      defaultBoundaryUID_;
  std::map< uint_t, BoundaryUID >  meshFlagToID_;
  std::map< BoundaryUID, DoFType > boundaryUIDToType_;

};



} // namespace hyteg