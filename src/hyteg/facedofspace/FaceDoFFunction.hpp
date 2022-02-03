/*
 * Copyright (c) 2017-2021 Daniel Drzisga, Dominik Thoennes, Marcus Mohr, Nils Kohl.
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

#include "hyteg/boundary/BoundaryConditions.hpp"
#include "hyteg/functions/Function.hpp"
#include "hyteg/functions/FunctionProperties.hpp"
#include "hyteg/p1functionspace/P1Function.hpp"

#include "FaceDoFDataHandling.hpp"
#include "FaceDoFMacroEdge.hpp"
#include "FaceDoFMacroFace.hpp"
#include "FaceDoFMacroVertex.hpp"
#include "FaceDoFMemory.hpp"
#include "FaceDoFPackInfo.hpp"

namespace hyteg {

template < typename ValueType >
class FaceDoFFunction : public Function< FaceDoFFunction< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = FaceDoFFunction< VType >;

   FaceDoFFunction( const std::string&                         name,
                    const std::shared_ptr< PrimitiveStorage >& storage,
                    uint_t                                     minLevel,
                    uint_t                                     maxLevel );

   FaceDoFFunction( const std::string&                         name,
                    const std::shared_ptr< PrimitiveStorage >& storage,
                    uint_t                                     minLevel,
                    uint_t                                     maxLevel,
                    BoundaryCondition                          boundaryCondition );

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType flag = All ) const;

   void interpolate( ValueType constant, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                     const std::vector< std::reference_wrapper< const FaceDoFFunction< ValueType > > >&   srcFunctions,
                     uint_t                                                                               level,
                     DoFType                                                                              flag = All ) const;

   void interpolate( const std::vector< std::function< ValueType( const hyteg::Point3D& ) > >& expr,
                     uint_t                                                                    level,
                     DoFType                                                                   flag = All ) const;

   void assign( const std::vector< ValueType >&                                                    scalars,
                const std::vector< std::reference_wrapper< const FaceDoFFunction< ValueType > > >& functions,
                uint_t                                                                             level,
                DoFType                                                                            flag = All ) const;

   void add( ValueType scalar, uint_t level, DoFType flag = All ) const;

   void add( const std::vector< ValueType >&                                                    scalars,
             const std::vector< std::reference_wrapper< const FaceDoFFunction< ValueType > > >& functions,
             uint_t                                                                             level,
             DoFType                                                                            flag = All ) const;

   void enumerate( uint_t level, ValueType offset );

   void enumerate( uint_t level );

   ValueType getMaxValue( uint_t level, DoFType flag = All );
   ValueType getMinValue( uint_t level, DoFType flag = All );
   ValueType getMaxMagnitude( uint_t level, DoFType flag = All );

   const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& getVertexDataID() const { return vertexDataID_; }

   const PrimitiveDataID< FunctionMemory< ValueType >, Edge >& getEdgeDataID() const { return edgeDataID_; }

   const PrimitiveDataID< FunctionMemory< ValueType >, Face >& getFaceDataID() const { return faceDataID_; }

   BoundaryCondition getBoundaryCondition() const { return boundaryCondition_; }

   template < typename SenderType, typename ReceiverType >
   void startCommunication( const uint_t& level ) const
   {
      communicators_.at( level )->template startCommunication< SenderType, ReceiverType >();
   }

   template < typename SenderType, typename ReceiverType >
   void endCommunication( const uint_t& level ) const
   {
      communicators_.at( level )->template endCommunication< SenderType, ReceiverType >();
   }

   template < typename SenderType, typename ReceiverType >
   void communicate( const uint_t& level ) const
   {
      startCommunication< SenderType, ReceiverType >( level );
      endCommunication< SenderType, ReceiverType >( level );
   }

   void setLocalCommunicationMode( const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode )
   {
      for ( auto& communicator : communicators_ )
      {
         communicator.second->setLocalCommunicationMode( localCommunicationMode );
      }
   }

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
   void multElementwise( const std::vector< std::reference_wrapper< const FaceDoFFunction< ValueType > > >& functions,
                         uint_t                                                                             level,
                         DoFType                                                                            flag = All ) const;

   void copyFrom( const FaceDoFFunction< ValueType >& other, const uint_t& level ) const;

   void copyFrom( const FaceDoFFunction< ValueType >&            other,
                  const uint_t&                                  level,
                  const std::map< PrimitiveID::IDType, uint_t >& localPrimitiveIDsToRank,
                  const std::map< PrimitiveID::IDType, uint_t >& otherPrimitiveIDsToRank ) const;

   void setBoundaryCondition( BoundaryCondition bc ) { boundaryCondition_ = std::move( bc ); }

   ValueType dotLocal( const FaceDoFFunction< ValueType >& secondOp, uint_t level, DoFType flag ) const;

   ValueType dotGlobal( const FaceDoFFunction< ValueType >& secondOp, uint_t level, DoFType flag ) const;

   void swap( const FaceDoFFunction< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const;

   /// conversion to/from linear algebra representation
   /// @{
   void toVector( const FaceDoFFunction< idx_t >&       numerator,
                  const std::shared_ptr< VectorProxy >& vec,
                  uint_t                                level,
                  DoFType                               flag ) const
   {
      WALBERLA_ABORT( "Congrats :( You have detected another unimplemented feature of FaceDoFFunction" );
   }

   void fromVector( const FaceDoFFunction< idx_t >&       numerator,
                    const std::shared_ptr< VectorProxy >& vec,
                    uint_t                                level,
                    DoFType                               flag ) const
   {
      WALBERLA_ABORT( "Congrats :( You have detected another unimplemented feature of FaceDoFFunction" );
   }
   /// @}

 private:
   using Function< FaceDoFFunction< ValueType > >::communicators_;

   PrimitiveDataID< FunctionMemory< ValueType >, Vertex > vertexDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Edge >   edgeDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Face >   faceDataID_;

   BoundaryCondition boundaryCondition_;
};

} // namespace hyteg
