/*
 * Copyright (c) 2017-2019 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "core/DataTypes.h"

#include "hyteg/Function.hpp"
#include "hyteg/edgedofspace/EdgeDoFFunction.hpp"
#include "hyteg/p1functionspace/VertexDoFFunction.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"
// #include "hyteg/p2functionspace/P2Multigrid.hpp"
// #include "hyteg/p2functionspace/P2TransferOperators.hpp"
// #include "hyteg/p2functionspace/P2MacroFace.hpp"

namespace hyteg {

using walberla::real_c;

template < typename ValueType >
class P2Function : public Function< P2Function< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = P2Function< VType >;

   P2Function( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage );

   P2Function( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage, uint_t minLevel, uint_t maxLevel );

   P2Function( const std::string&                         name,
               const std::shared_ptr< PrimitiveStorage >& storage,
               uint_t                                     minLevel,
               uint_t                                     maxLevel,
               BoundaryCondition                          boundaryCondition );

   inline vertexdof::VertexDoFFunction< ValueType > getVertexDoFFunctionCopy() const { return vertexDoFFunction_; }
   inline EdgeDoFFunction< ValueType >              getEdgeDoFFunctionCopy() const { return edgeDoFFunction_; }

   inline const vertexdof::VertexDoFFunction< ValueType >& getVertexDoFFunction() const { return vertexDoFFunction_; }
   inline const EdgeDoFFunction< ValueType >&              getEdgeDoFFunction() const { return edgeDoFFunction_; }

   template < typename SenderType, typename ReceiverType >
   void communicate( const uint_t& level ) const
   {
      vertexDoFFunction_.template communicate< SenderType, ReceiverType >( level );
      edgeDoFFunction_.template communicate< SenderType, ReceiverType >( level );
   }

   inline ValueType evaluate( const Point3D& coordinates, uint_t level ) const;

   inline void evaluateGradient( const Point3D& coordinates, uint_t level, Point3D& gradient ) const;

   void interpolate( const ValueType& constant, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, BoundaryUID boundaryUID ) const;

   void interpolate( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                     const std::vector< std::reference_wrapper< const P2Function< ValueType > > >&        srcFunctions,
                     uint_t                                                                               level,
                     DoFType                                                                              flag = All ) const;

   void swap( const P2Function< ValueType >& other, const uint_t& level, const DoFType& dofType = All ) const;

   void assign( const std::vector< ValueType >&                                               scalars,
                const std::vector< std::reference_wrapper< const P2Function< ValueType > > >& functions,
                uint_t                                                                        level,
                DoFType                                                                       flag = All ) const;

   void assign( const P1Function< ValueType >& src, const uint_t& P2Level, const DoFType& flag = All ) const;

   void add( const ValueType& scalar, uint_t level, DoFType flag = All ) const;

   void add( const std::vector< ValueType >&                                               scalars,
             const std::vector< std::reference_wrapper< const P2Function< ValueType > > >& functions,
             uint_t                                                                        level,
             DoFType                                                                       flag = All ) const;

   ValueType dotGlobal( const P2Function< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const;

   ValueType dotLocal( const P2Function< ValueType >& rhs, uint_t level, const DoFType& flag = All ) const;

   ValueType sumGlobal( uint_t level, const DoFType& flag = All, const bool & absolute = false ) const;

   ValueType sumLocal( uint_t level, const DoFType& flag = All, const bool & absolute = false ) const;

   void prolongateP1ToP2( const hyteg::P1Function< ValueType >& p1Function, const uint_t& level, const DoFType& flag = All ) const;

   void restrictP2ToP1( const P1Function< ValueType >& p1Function, const uint_t& level, const DoFType& flag = All ) const;

   void restrictInjection( uint_t sourceLevel, DoFType flag = All ) const;

   ValueType getMaxValue( uint_t level, DoFType flag = All ) const;

   ValueType getMaxMagnitude( uint_t level, DoFType flag = All ) const;

   ValueType getMinValue( uint_t level, DoFType flag = All ) const;

   BoundaryCondition getBoundaryCondition() const;

   void enumerate( uint_t level ) const;

   void enumerate( uint_t level, ValueType& offset ) const;

   void setLocalCommunicationMode( const communication::BufferedCommunicator::LocalCommunicationMode& localCommMode );

 private:
   using Function< P2Function< ValueType > >::communicators_;

   vertexdof::VertexDoFFunction< ValueType > vertexDoFFunction_;
   EdgeDoFFunction< ValueType >              edgeDoFFunction_;
};

namespace p2function {

void projectMean( const P2Function< real_t >& pressure, const uint_t& level );

} // namespace p2function

} //namespace hyteg
