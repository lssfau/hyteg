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

#include "core/mpi/all.h"

#include "hyteg/Function.hpp"
#include "hyteg/boundary/BoundaryConditions.hpp"

namespace hyteg {

using walberla::real_c;

template< typename ValueType >
class FunctionMemory;
template< typename DataType, typename PrimitiveType >
class PrimitiveDataID;

class PrimitiveStorage;
class Vertex;
class Edge;
class Face;
class Cell;

namespace edgedof {
///@name Size Functions
///@{

uint_t edgeDoFMacroVertexFunctionMemorySize( const uint_t &level, const Primitive & primitive );

uint_t edgeDoFMacroEdgeFunctionMemorySize( const uint_t &level, const Primitive & primitive );

uint_t edgeDoFMacroFaceFunctionMemorySize( const uint_t &level, const Primitive & primitive );

uint_t edgeDoFMacroCellFunctionMemorySize( const uint_t & level, const Primitive & primitive );

///@}

}// namespace edgedof

template< typename ValueType >
class EdgeDoFFunction : public Function< EdgeDoFFunction< ValueType > >
{
public:
  EdgeDoFFunction( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage );

  EdgeDoFFunction( const std::string&                         name,
                   const std::shared_ptr< PrimitiveStorage >& storage,
                   const uint_t&                              minLevel,
                   const uint_t&                              maxLevel );

  EdgeDoFFunction( const std::string&                         name,
                   const std::shared_ptr< PrimitiveStorage >& storage,
                   const uint_t&                              minLevel,
                   const uint_t&                              maxLevel,
                   const BoundaryCondition&                   boundaryCondition,
                   const DoFType&                             boundaryTypeToSkipDuringAdditiveCommunication = DoFType::DirichletBoundary );

  void swap( const EdgeDoFFunction< ValueType > & other,
             const uint_t & level,
             const DoFType & flag = All ) const;

  void assign( const std::vector< ValueType >&                                                    scalars,
               const std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >& functions,
               uint_t                                                                             level,
               DoFType                                                                            flag = All ) const;

  void add( const ValueType& scalar, uint_t level, DoFType flag = All ) const;

  void add( const std::vector< ValueType >&                                                    scalars,
            const std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >& functions,
            uint_t  level,
            DoFType flag = All ) const;

  /// Interpolates a given expression to a EdgeDoFFunction

  void interpolate( const ValueType& constant, uint_t level, DoFType flag = All ) const;

  void interpolate( const std::function< ValueType( const Point3D & ) >& expr,
                          uint_t level, DoFType flag = All) const;

  void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, BoundaryUID boundaryUID ) const;

  void interpolateExtended( const std::function<ValueType(const Point3D &, const std::vector<ValueType>&)> &expr,
                            const std::vector<EdgeDoFFunction<ValueType>*> srcFunctions,
                            uint_t level,
                            DoFType flag = All) const;

  void interpolateExtended( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                            const std::vector< EdgeDoFFunction< ValueType >* >                                   srcFunctions,
                            uint_t                                                                               level,
                            BoundaryUID                                                                          boundaryUID ) const;

   /// Compute the product of several functions in an elementwise fashion
   ///
   /// The method takes as input a collection of functions. These are multiplied together in an elementwise fashion.
   /// The latter is to be understood not in a FE context, but in the sense of element-wise operators in matrix/array
   /// oriented languages, i.e. the product is a function of the same type as the inputs and its DoFs are formed as
   /// product of the corresponding DoFs of the input functions. The result is stored in the function object on which
   /// the method is invoked, overwritting its contents. It is safe, if the destination function is part of the product.
   ///
   /// \param functions  the functions forming the product
   /// \param level      level on which the multiplication should be computed
   /// \param flag       marks those primitives which are partaking in the computation of the product
   void multElementwise( const std::vector< std::reference_wrapper< const EdgeDoFFunction< ValueType > > >& functions,
                         uint_t                                                                               level,
                         DoFType                                                                              flag = All ) const;

  ValueType dotLocal(const EdgeDoFFunction <ValueType> &rhs, const uint_t level, const DoFType flag = All) const;

  ValueType sumLocal( const uint_t& level, const DoFType& flag = All, const bool & absolute = false ) const;
  ValueType sumGlobal( const uint_t& level, const DoFType& flag = All, const bool & absolute = false ) const;

  void enumerate( uint_t level ) const;

  const PrimitiveDataID< FunctionMemory< ValueType >, Vertex>   & getVertexDataID() const { return vertexDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >,   Edge>   & getEdgeDataID()   const { return edgeDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >,   Face>   & getFaceDataID()   const { return faceDataID_; }
  const PrimitiveDataID< FunctionMemory< ValueType >,   Cell>   & getCellDataID()   const { return cellDataID_; }


  ValueType getMaxValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;
  ValueType getMinValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;
  ValueType getMaxMagnitude( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;

  inline BoundaryCondition getBoundaryCondition() const { return boundaryCondition_; }
  inline DoFType           getBoundaryTypeToSkipDuringAdditiveCommunication() const { return boundaryTypeToSkipDuringAdditiveCommunication_; }

  template< typename SenderType, typename ReceiverType >
  inline void startCommunication( const uint_t & level ) const
  {
    if ( isDummy() ) { return; }
    communicators_.at( level )->template startCommunication< SenderType, ReceiverType >();
  }

  template< typename SenderType, typename ReceiverType >
  inline void endCommunication( const uint_t & level ) const
  {
    if ( isDummy() ) { return; }
    communicators_.at( level )->template endCommunication< SenderType, ReceiverType >();
  }

  template< typename SenderType, typename ReceiverType >
  inline void communicate( const uint_t & level ) const
  {
    if ( isDummy() ) { return; }
    communicators_.at( level )->template communicate< SenderType, ReceiverType >();
  }

  template < typename SenderType, typename ReceiverType >
  inline void startAdditiveCommunication( const uint_t& level ) const
  {
    if( isDummy() )
    {
      return;
    }
    interpolateByPrimitiveType< ReceiverType >(
    real_c( 0 ), level, DoFType::All ^ boundaryTypeToSkipDuringAdditiveCommunication_ );
    additiveCommunicators_.at( level )->template startCommunication< SenderType, ReceiverType >();
  }

  template < typename SenderType, typename ReceiverType >
  inline void endAdditiveCommunication( const uint_t& level ) const
  {
    if( isDummy() )
    {
      return;
    }
    additiveCommunicators_.at( level )->template endCommunication< SenderType, ReceiverType >();
  }

  template < typename SenderType, typename ReceiverType >
  inline void communicateAdditively( const uint_t& level ) const
  {
    if( isDummy() )
    {
      return;
    }
    interpolateByPrimitiveType< ReceiverType >(
    real_c( 0 ), level, DoFType::All ^ boundaryTypeToSkipDuringAdditiveCommunication_ );
    additiveCommunicators_.at( level )->template communicate< SenderType, ReceiverType >();
  }

  inline void setLocalCommunicationMode( const communication::BufferedCommunicator::LocalCommunicationMode & localCommunicationMode )
  {
    if ( isDummy() ) { return; }
    for ( auto & communicator : communicators_ )
    {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
    }
    for( auto & communicator : additiveCommunicators_ )
    {
      communicator.second->setLocalCommunicationMode( localCommunicationMode );
    }
  }

   using Function< EdgeDoFFunction< ValueType > >::isDummy;

private:

   template < typename PrimitiveType >
   void interpolateByPrimitiveType( const ValueType& constant, uint_t level, DoFType flag = All ) const;

   void enumerate( uint_t level, ValueType& offset ) const;

   using Function< EdgeDoFFunction< ValueType > >::communicators_;
   using Function< EdgeDoFFunction< ValueType > >::additiveCommunicators_;

   PrimitiveDataID< FunctionMemory< ValueType >, Vertex > vertexDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Edge > edgeDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Face > faceDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Cell > cellDataID_;

   BoundaryCondition boundaryCondition_;

   DoFType boundaryTypeToSkipDuringAdditiveCommunication_;

   /// friend P2Function for usage of enumerate
   friend class P2Function< ValueType >;
};


}// namespace hyteg
