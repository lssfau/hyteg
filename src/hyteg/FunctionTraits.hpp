#pragma once

#include <string>

namespace hyteg {

// To define correct method signatures in the Function base class, we
// need to know the value type of the derived classes.
//
// (AFAIK) It is not possible to access typedefs of the derived class from the
// base class when using the CRTP. This means we need to provide them via traits to
// the function space base class "Function".
//
// See also: https://stackoverflow.com/a/6006629


/////////////////////
// Tag dispatching //
/////////////////////

struct VertexDoFFunctionTag {};
typedef VertexDoFFunctionTag P1FunctionTag;
struct EdgeDoFFunctionTag {};
struct DGFunctionTag {};
struct P2FunctionTag {};
struct P1StokesFunctionTag {};
struct P2P1TaylorHoodFunctionTag {};
struct P2P2StokesFunctionTag {};


//////////////////////////
// Forward declarations //
//////////////////////////

namespace vertexdof {

template< typename VType >
class VertexDoFFunction;

} // namespace vertexdof

template< typename VType >
class EdgeDoFFunction;

template< typename VType >
class DGFunction;

// Composites

template< typename VType >
class P2Function;

template< typename VType >
class P1StokesFunction;

template< typename VType >
class P2P1TaylorHoodFunction;

template< typename VType >
class P2P2StokesFunction;


///////////////////////////////////////////////////////////////////
// Function trait defining the value type of the derived classes //
///////////////////////////////////////////////////////////////////

/// Empty trait
template< typename FunctionType >
struct FunctionTrait;

/// Vertex DoF specialization
template< typename VType >
struct FunctionTrait< vertexdof::VertexDoFFunction< VType > >
{
  typedef VType ValueType;
  typedef VertexDoFFunctionTag Tag;

  static std::string getTypeName() { return "P1Function / VertexDoFFunction"; }
};

/// Edge DoF specialization
template< typename VType >
struct FunctionTrait< EdgeDoFFunction< VType > >
{
  typedef VType ValueType;
  typedef EdgeDoFFunctionTag Tag;

  static std::string getTypeName() { return "EdgeDoFFunction"; }
};

/// DG specialization
template< typename VType >
struct FunctionTrait< DGFunction< VType > >
{
  typedef VType ValueType;
  typedef DGFunctionTag Tag;

  static std::string getTypeName() { return "DGFunction"; }
};

/// P2 specialization
template< typename VType >
struct FunctionTrait< P2Function< VType > >
{
  typedef VType ValueType;
  typedef P2FunctionTag Tag;

  static std::string getTypeName() { return "P2Function"; }
};

/// P1Stokes specialization
template< typename VType >
struct FunctionTrait< P1StokesFunction< VType > >
{
    typedef VType ValueType;
    typedef P1StokesFunctionTag Tag;

    static std::string getTypeName() { return "P1StokesFunction"; }
};

/// P2P1TaylorHood specialization
template< typename VType >
struct FunctionTrait< P2P1TaylorHoodFunction< VType > >
{
    typedef VType ValueType;
    typedef P2P1TaylorHoodFunctionTag Tag;

    static std::string getTypeName() { return "P2P1TaylorHoodFunction"; }
};

/// P2P2Stokes specialization
template< typename VType >
struct FunctionTrait< P2P2StokesFunction< VType > >
{
    typedef VType ValueType;
    typedef P2P2StokesFunctionTag Tag;

    static std::string getTypeName() { return "P2P2StokesFunction"; }
};

}
