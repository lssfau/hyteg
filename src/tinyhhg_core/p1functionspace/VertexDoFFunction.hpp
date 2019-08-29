#pragma once

#include <memory>
#include <string>
#include <vector>

#include "core/DataTypes.h"

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/boundary/BoundaryConditions.hpp"
#include "tinyhhg_core/types/flags.hpp"
//TODO this should be improved but we need the enum which cant be forward declared
#include "tinyhhg_core/communication/BufferedCommunication.hpp"

namespace hhg {

using walberla::real_c;
using walberla::real_t;
using walberla::uint_t;

template < typename ValueType >
class DGFunction;
template < typename ValueType >
class FunctionMemory;
template < typename DataType, typename PrimitiveType >
class PrimitiveDataID;

class PrimitiveStorage;
class Vertex;
class Edge;
class Face;
class Cell;

namespace vertexdof {

template < typename ValueType >
class VertexDoFFunction : public Function< VertexDoFFunction< ValueType > >
{
 public:
   typedef ValueType valueType;

   template < typename VType >
   using FunctionType = vertexdof::VertexDoFFunction< VType >;

   VertexDoFFunction( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage );

   VertexDoFFunction( const std::string&                         name,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      uint_t                                     minLevel,
                      uint_t                                     maxLevel );

   VertexDoFFunction( const std::string&                         name,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      uint_t                                     minLevel,
                      uint_t                                     maxLevel,
                      BoundaryCondition                          boundaryCondition,
                      const DoFType& boundaryTypeToSkipDuringAdditiveCommunication = DoFType::DirichletBoundary );

   const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& getVertexDataID() const { return vertexDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&   getEdgeDataID() const { return edgeDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Face >&   getFaceDataID() const { return faceDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&   getCellDataID() const { return cellDataID_; }

   void swap( const VertexDoFFunction< ValueType >& other, const uint_t& level, const DoFType& flag = All ) const;

   real_t evaluate( const Point3D& coordinates, uint_t level ) const;

   void evaluateGradient( const Point3D& coordinates, uint_t level, Point3D& gradient ) const;

   void assign( const std::vector< ValueType >&                                                      scalars,
                const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
                uint_t                                                                               level,
                DoFType                                                                              flag = All ) const;

   void assign( const P2Function< ValueType >& src, const uint_t& P1Level, const DoFType& flag = All ) const;

   void add( const ValueType& scalar, const uint_t& level, DoFType flag = All ) const;

   void add( const std::vector< ValueType >&                                                      scalars,
             const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
             uint_t                                                                               level,
             DoFType                                                                              flag = All ) const;

   void multElementwise( const std::vector< std::reference_wrapper< const VertexDoFFunction< ValueType > > >& functions,
                         uint_t                                                                               level,
                         DoFType                                                                              flag = All ) const;

   ValueType dotLocal( const VertexDoFFunction< ValueType >& rhs, uint_t level, DoFType flag = All ) const;
   ValueType dotGlobal( const VertexDoFFunction< ValueType >& rhs, uint_t level, DoFType flag = All ) const;

   ValueType sumLocal( const uint_t& level, const DoFType& flag = All, const bool & absolute = false ) const;
   ValueType sumGlobal( const uint_t& level, const DoFType& flag = All, const bool & absolute = false ) const;

   void integrateDG( DGFunction< ValueType >& rhs, VertexDoFFunction< ValueType >& rhsP1, uint_t level, DoFType flag );

   /// Interpolates a given expression to a VertexDoFFunction
   void interpolate( const ValueType& constant, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, BoundaryUID boundaryUID ) const;

   void interpolateExtended( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                             const std::vector< VertexDoFFunction* >                                              srcFunctions,
                             uint_t                                                                               level,
                             DoFType flag = All ) const;

   void interpolateExtended( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                             const std::vector< VertexDoFFunction* >                                              srcFunctions,
                             uint_t                                                                               level,
                             BoundaryUID boundaryUID ) const;

   /// assigns unique values to all data points
   /// this function is mainly used for petsc to get global identifier for all DoFs
   /// \tparam ValueType
   /// \param level
   void enumerate( uint_t level ) const;

   // TODO: write more general version(s)
   ValueType getMaxValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;
   ValueType getMinValue( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;
   ValueType getMaxMagnitude( uint_t level, DoFType flag = All, bool mpiReduce = true ) const;

   BoundaryCondition getBoundaryCondition() const;
   inline DoFType    getBoundaryTypeToSkipDuringAdditiveCommunication() const
   {
      return boundaryTypeToSkipDuringAdditiveCommunication_;
   }

   template < typename SenderType, typename ReceiverType >
   inline void startCommunication( const uint_t& level ) const
   {
      if ( isDummy() )
      {
         return;
      }
      communicators_.at( level )->template startCommunication< SenderType, ReceiverType >();
   }

   template < typename SenderType, typename ReceiverType >
   inline void endCommunication( const uint_t& level ) const
   {
      if ( isDummy() )
      {
         return;
      }
      communicators_.at( level )->template endCommunication< SenderType, ReceiverType >();
   }

   template < typename SenderType, typename ReceiverType >
   inline void communicate( const uint_t& level ) const
   {
      if ( isDummy() )
      {
         return;
      }
      communicators_.at( level )->template communicate< SenderType, ReceiverType >();
   }

   template < typename SenderType, typename ReceiverType >
   inline void startAdditiveCommunication( const uint_t& level ) const
   {
      if ( isDummy() )
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
      if ( isDummy() )
      {
         return;
      }
      additiveCommunicators_.at( level )->template endCommunication< SenderType, ReceiverType >();
   }

   template < typename SenderType, typename ReceiverType >
   inline void communicateAdditively( const uint_t& level ) const
   {
      if ( isDummy() )
      {
         return;
      }
      interpolateByPrimitiveType< ReceiverType >(
          real_c( 0 ), level, DoFType::All ^ boundaryTypeToSkipDuringAdditiveCommunication_ );
      additiveCommunicators_.at( level )->template communicate< SenderType, ReceiverType >();
   }

   void setLocalCommunicationMode( const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode );

   using Function< VertexDoFFunction< ValueType > >::isDummy;

 private:
   template < typename PrimitiveType >
   void interpolateByPrimitiveType( const ValueType& constant, uint_t level, DoFType flag = All ) const;

   void enumerate( uint_t level, ValueType& offset ) const;

   using Function< VertexDoFFunction< ValueType > >::communicators_;
   using Function< VertexDoFFunction< ValueType > >::additiveCommunicators_;

   PrimitiveDataID< FunctionMemory< ValueType >, Vertex > vertexDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Edge >   edgeDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Face >   faceDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Cell >   cellDataID_;

   BoundaryCondition boundaryCondition_;

   DoFType boundaryTypeToSkipDuringAdditiveCommunication_;

   /// friend Stokes and P2Function for usage of enumerate
   friend class P2Function< ValueType >;
   friend class P1StokesFunction< ValueType >;
   friend class P2P1TaylorHoodFunction< ValueType >;
};

inline void projectMean( const VertexDoFFunction< real_t >& pressure, const uint_t& level )
{
   if ( pressure.isDummy() )
   {
      return;
   }
   const uint_t numGlobalVertices = numberOfGlobalDoFs< VertexDoFFunctionTag >( *pressure.getStorage(), level );
   const real_t sum               = pressure.sumGlobal( level, All );
   pressure.add( -sum / real_c( numGlobalVertices ), level, All );
}

} // namespace vertexdof
} // namespace hhg
