#pragma once

#include "core/DataTypes.h"

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/edgedofspace/EdgeDoFFunction.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFFunction.hpp"
#include "tinyhhg_core/p1functionspace/P1Function.hpp"

namespace hhg {

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
               BoundaryCondition                          boundaryCondition,
               const DoFType& boundaryTypeToSkipDuringAdditiveCommunication = DoFType::DirichletBoundary );

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

   void interpolate( const ValueType& constant, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, BoundaryUID boundaryUID ) const;

   void interpolateExtended( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                             const std::vector< P2Function< ValueType >* >                                        srcFunctions,
                             uint_t                                                                               level,
                             DoFType flag = All ) const;

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

   void prolongateP1ToP2( const hhg::P1Function< ValueType >& p1Function, const uint_t& level, const DoFType& flag = All ) const;

   void restrictP2ToP1( const P1Function< ValueType >& p1Function, const uint_t& level, const DoFType& flag = All ) const;

   void restrictInjection( uint_t sourceLevel, DoFType flag = All ) const;

   void interpolate( std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                     const std::vector< P2Function< ValueType >* >                                  srcFunctions,
                     uint_t                                                                         level,
                     DoFType                                                                        flag = All ) const;

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

} //namespace hhg
