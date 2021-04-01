/*
 * Copyright (c) 2017-2019 Dominik Thoennes, Nils Kohl.
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
struct P1VectorFunctionTag {};
struct P1VectorFunctionTag_AltKind {};
struct P2VectorFunctionTag {};


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

template< typename VType >
class P1VectorFunction;

template< typename VType >
class P1VectorFunction_AltKind;

template< typename VType >
class P2VectorFunction;


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

/// P1VectorFunction specialization
template< typename VType >
struct FunctionTrait< P1VectorFunction< VType > >
{
   typedef VType ValueType;
   typedef P1VectorFunctionTag Tag;

   static std::string getTypeName() { return "P1VectorFunction"; }
};

/// P1VectorFunction specialization (AltKind)
template< typename VType >
struct FunctionTrait< P1VectorFunction_AltKind< VType > >
{
   typedef VType ValueType;
   typedef P1VectorFunctionTag_AltKind Tag;

   static std::string getTypeName() { return "P1VectorFunction_AltKind"; }
};

/// P2VectorFunction specialization
template< typename VType >
struct FunctionTrait< P2VectorFunction< VType > >
{
   typedef VType ValueType;
   typedef P2VectorFunctionTag Tag;

   static std::string getTypeName() { return "P2VectorFunction"; }
};

}
