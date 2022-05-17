/*
 * Copyright (c) 2017-2021 Dominik Thoennes, Nils Kohl, Marcus Mohr.
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

struct P0FunctionTag
{};
struct VertexDoFFunctionTag
{};
typedef VertexDoFFunctionTag P1FunctionTag;
struct EdgeDoFFunctionTag
{};
struct VolumeFunctionTag
{};
struct FaceDoFFunction_old_Tag
{};
struct DGFunctionTag
{};
struct P2FunctionTag
{};
struct P1StokesFunctionTag
{};
struct P2P1TaylorHoodFunctionTag
{};
struct P2P1TaylorHoodBlockFunctionTag
{};
struct P2P2StokesFunctionTag
{};
struct P1VectorFunctionTag
{};
struct P2VectorFunctionTag
{};
struct DGVectorFunctionTag
{};
struct P1DGEFunctionTag
{};
struct P1DGEP1StokesFunctionTag
{};

//////////////////////////
// Forward declarations //
//////////////////////////

template < typename VType >
class P0Function;

namespace vertexdof {

template < typename VType >
class VertexDoFFunction;

} // namespace vertexdof

template < typename VType >
class EdgeDoFFunction;

namespace volumedofspace {
template < typename VType >
class VolumeDoFFunction;
}

template < typename VType >
class FaceDoFFunction_old;

namespace dg {
template < typename VType >
class DGFunction;
}

// Composites

template < typename VType >
class P2Function;

template < typename VType >
class P1StokesFunction;

template < typename VType >
class P2P1TaylorHoodFunction;

template < typename VType >
class P2P1TaylorHoodBlockFunction;

template < typename VType >
class P2P2StokesFunction;

template < typename VType >
class P1VectorFunction;

template < typename VType >
class P1VectorFunction_AltKind;

template < typename VType >
class P2VectorFunction;

namespace dg {
template < typename VType >
class DGVectorFunction;
}

template < typename VType >
class P1DGEFunction;

template < typename VType >
class P1DGEP0StokesFunction;

///////////////////////////////////////////////////////////////////
// Enum for getting info on type of a GenericFunction
///////////////////////////////////////////////////////////////////

namespace functionTraits {

typedef enum
{
   P0_FUNCTION,
   P1_FUNCTION,
   P2_FUNCTION,
   EDGE_DOF_FUNCTION,
   VOLUME_DOF_FUNCTION,
   FACE_DOF_FUNCTION_OLD,
   DG_FUNCTION,
   P1_VECTOR_FUNCTION,
   P2_VECTOR_FUNCTION,
   DG_VECTOR_FUNCTION,
   P1DGE_FUNCTION,
   OTHER_FUNCTION
} FunctionKind;

}

///////////////////////////////////////////////////////////////////
// Function trait defining the value type of the derived classes //
///////////////////////////////////////////////////////////////////

/// Empty trait
template < typename FunctionType >
struct FunctionTrait;

/// P0 DoF specialization
template < typename VType >
struct FunctionTrait< P0Function< VType > >
{
   typedef VType         ValueType;
   typedef P0FunctionTag Tag;

   static std::string getTypeName() { return "P0Function"; }

   static const functionTraits::FunctionKind kind = functionTraits::P0_FUNCTION;
};

/// Vertex DoF specialization
template < typename VType >
struct FunctionTrait< vertexdof::VertexDoFFunction< VType > >
{
   typedef VType                     ValueType;
   typedef VertexDoFFunctionTag      Tag;
   typedef P1VectorFunction< VType > AssocVectorFunctionType;

   static std::string getTypeName() { return "P1Function / VertexDoFFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::P1_FUNCTION;
};

/// Edge DoF specialization
template < typename VType >
struct FunctionTrait< EdgeDoFFunction< VType > >
{
   typedef VType              ValueType;
   typedef EdgeDoFFunctionTag Tag;

   static std::string getTypeName() { return "EdgeDoFFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::EDGE_DOF_FUNCTION;
};

/// Face DoF specialization
template < typename VType >
struct FunctionTrait< FaceDoFFunction_old< VType > >
{
   typedef VType                   ValueType;
   typedef FaceDoFFunction_old_Tag Tag;

   static std::string getTypeName() { return "FaceDoFFunction_old"; }

   static const functionTraits::FunctionKind kind = functionTraits::FACE_DOF_FUNCTION_OLD;
};

/// Volume DoF specialization
template < typename VType >
struct FunctionTrait< volumedofspace::VolumeDoFFunction< VType > >
{
   typedef VType         ValueType;
   typedef DGFunctionTag Tag;

   static std::string getTypeName() { return "VolumeDoFFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::VOLUME_DOF_FUNCTION;
};

/// DG specialization
template < typename VType >
struct FunctionTrait< dg::DGFunction< VType > >
{
   typedef VType         ValueType;
   typedef DGFunctionTag Tag;

   static std::string getTypeName() { return "DGFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::DG_FUNCTION;
};

/// P2 specialization
template < typename VType >
struct FunctionTrait< P2Function< VType > >
{
   typedef VType                     ValueType;
   typedef P2FunctionTag             Tag;
   typedef P2VectorFunction< VType > AssocVectorFunctionType;

   static std::string getTypeName() { return "P2Function"; }

   static const functionTraits::FunctionKind kind = functionTraits::P2_FUNCTION;
};

/// P1Stokes specialization
template < typename VType >
struct FunctionTrait< P1StokesFunction< VType > >
{
   typedef VType               ValueType;
   typedef P1StokesFunctionTag Tag;

   static std::string getTypeName() { return "P1StokesFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::OTHER_FUNCTION;
};

/// P2P1TaylorHood specialization
template < typename VType >
struct FunctionTrait< P2P1TaylorHoodFunction< VType > >
{
   typedef VType                     ValueType;
   typedef P2P1TaylorHoodFunctionTag Tag;

   static std::string getTypeName() { return "P2P1TaylorHoodFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::OTHER_FUNCTION;
};

/// P2P1TaylorHood specialization
template < typename VType >
struct FunctionTrait< P2P1TaylorHoodBlockFunction< VType > >
{
   typedef VType                          ValueType;
   typedef P2P1TaylorHoodBlockFunctionTag Tag;

   static std::string getTypeName() { return "P2P1TaylorHoodBlockFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::OTHER_FUNCTION;
};

/// P2P2Stokes specialization
template < typename VType >
struct FunctionTrait< P2P2StokesFunction< VType > >
{
   typedef VType                 ValueType;
   typedef P2P2StokesFunctionTag Tag;

   static std::string getTypeName() { return "P2P2StokesFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::OTHER_FUNCTION;
};

/// P1VectorFunction specialization
template < typename VType >
struct FunctionTrait< P1VectorFunction< VType > >
{
   typedef VType                                 ValueType;
   typedef P1VectorFunctionTag                   Tag;
   typedef vertexdof::VertexDoFFunction< VType > VectorComponentType;
   // I'd prefer to have P1Function< VType > above, but couldn't get
   // that to work. Forward declaration of P1Function as class fails
   // since it is no class, see P1Function.hpp

   static std::string getTypeName() { return "P1VectorFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::P1_VECTOR_FUNCTION;
};

/// P2VectorFunction specialization
template < typename VType >
struct FunctionTrait< P2VectorFunction< VType > >
{
   typedef VType               ValueType;
   typedef P2VectorFunctionTag Tag;
   typedef P2Function< VType > VectorComponentType;

   static std::string getTypeName() { return "P2VectorFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::P2_VECTOR_FUNCTION;
};

/// DGVectorFunction specialization
template < typename VType >
struct FunctionTrait< dg::DGVectorFunction< VType > >
{
   typedef VType                   ValueType;
   typedef DGVectorFunctionTag     Tag;
   typedef dg::DGFunction< VType > VectorComponentType;

   static std::string getTypeName() { return "DGVectorFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::DG_VECTOR_FUNCTION;
};

/// P1DGEFunction specialization
template < typename VType >
struct FunctionTrait< P1DGEFunction< VType > >
{
   typedef VType                  ValueType;
   typedef P1DGEFunctionTag       Tag;
   typedef P1DGEFunction< VType > VectorComponentType;

   static std::string getTypeName() { return "P1DGEVectorFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::P1DGE_FUNCTION;
};

/// P1DGEP1StokesFunction specialization
template < typename VType >
struct FunctionTrait< P1DGEP0StokesFunction< VType > >
{
   typedef VType                    ValueType;
   typedef P1DGEP1StokesFunctionTag Tag;

   static std::string getTypeName() { return "P1DGEP1StokesFunction"; }

   static const functionTraits::FunctionKind kind = functionTraits::OTHER_FUNCTION;
};

} // namespace hyteg
