
#pragma once

namespace hhg {

// To define correct method signatures in the Function base class, we
// need to know the value type of the derived classes.
//
// (AFAIK) It is not possible to access typedefs of the derived class from the
// base class when using the CRTP. This means we need to provide them via traits to
// the function space base class "Function".
//
// See also: https://stackoverflow.com/a/6006629

//////////////////////////
// Forward declarations //
//////////////////////////

template< typename VType >
class P1Function;

template< typename VType >
class EdgeDoFFunction;

template< typename VType >
class BubbleFunction;

template< typename VType >
class DGFunction;


///////////////////////////////////////////////////////////////////
// Function trait defining the value type of the derived classes //
///////////////////////////////////////////////////////////////////

// Empty trait

template< typename FunctionType >
struct FunctionTrait;

// P1 specialization (Vertex DoFs)

template< typename VType >
struct FunctionTrait< P1Function< VType > >
{
  typedef VType ValueType;

  static std::string getTypeName() { return "P1Function"; }
};

// Edge Dof specialization
template< typename VType >
struct FunctionTrait< EdgeDoFFunction< VType > >
{
  typedef VType ValueType;

  static std::string getTypeName() { return "EdgeDoFFunction"; }
};

// Bubble specialization (Cell DoFs)

template< typename VType >
struct FunctionTrait< BubbleFunction< VType > >
{
  typedef VType ValueType;

  static std::string getTypeName() { return "BubbleFunction"; }
};

// DG specialization

template< typename VType >
struct FunctionTrait< DGFunction< VType > >
{
  typedef VType ValueType;

  static std::string getTypeName() { return "DGFunction"; }
};


}
