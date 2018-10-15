#pragma once

#include "tinyhhg_core/Function.hpp"
#include "tinyhhg_core/FunctionMemory.hpp"
#include "tinyhhg_core/FunctionProperties.hpp"
#include "tinyhhg_core/PrimitiveID.hpp"
#include "tinyhhg_core/boundary/BoundaryConditions.hpp"
#include "tinyhhg_core/communication/BufferedCommunication.hpp"
#include "tinyhhg_core/p1functionspace/P1DataHandling.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFAdditivePackInfo.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroCell.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroEdge.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroFace.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMacroVertex.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFMemory.hpp"
#include "tinyhhg_core/p1functionspace/VertexDoFPackInfo.hpp"
#include "tinyhhg_core/primitivedata/PrimitiveDataID.hpp"
#include "tinyhhg_core/primitives/Edge.hpp"
#include "tinyhhg_core/primitives/Face.hpp"
#include "tinyhhg_core/primitives/Vertex.hpp"
#include "tinyhhg_core/types/pointnd.hpp"

namespace hhg {
namespace vertexdof {

template < typename ValueType >
class VertexDoFFunction : public Function< VertexDoFFunction< ValueType > >
{
 public:
   VertexDoFFunction( const std::string& name, const std::shared_ptr< PrimitiveStorage >& storage );

   VertexDoFFunction( const std::string&                         name,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      uint_t                                     minLevel,
                      uint_t                                     maxLevel );

   VertexDoFFunction( const std::string&                         name,
                      const std::shared_ptr< PrimitiveStorage >& storage,
                      uint_t                                     minLevel,
                      uint_t                                     maxLevel,
                      BoundaryCondition                          boundaryCondition );
   const PrimitiveDataID< FunctionMemory< ValueType >, Vertex >& getVertexDataID() const { return vertexDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Edge >&   getEdgeDataID() const { return edgeDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Face >&   getFaceDataID() const { return faceDataID_; }
   const PrimitiveDataID< FunctionMemory< ValueType >, Cell >&   getCellDataID() const { return cellDataID_; }

   void assign( const std::vector< ValueType >                       scalars,
                const std::vector< VertexDoFFunction< ValueType >* > functions,
                uint_t                                               level,
                DoFType                                              flag = All );

   void add( const ValueType& scalar, const uint_t& level, DoFType flag = All );

   void add( const std::vector< ValueType >                       scalars,
             const std::vector< VertexDoFFunction< ValueType >* > functions,
             uint_t                                               level,
             DoFType                                              flag = All );

   void multElementwise( const std::vector< VertexDoFFunction< ValueType >* > functions, uint_t level, DoFType flag = All );

   real_t dotLocal( VertexDoFFunction< ValueType >& rhs, uint_t level, DoFType flag = All );
   real_t dotGlobal( VertexDoFFunction< ValueType >& rhs, uint_t level, DoFType flag = All );

   void integrateDG( DGFunction< ValueType >& rhs, VertexDoFFunction< ValueType >& rhsP1, uint_t level, DoFType flag );

   /// Interpolates a given expression to a VertexDoFFunction
   void interpolate( const ValueType& constant, uint_t level, DoFType flag = All ) const;

   void interpolate( const std::function< ValueType( const Point3D& ) >& expr, uint_t level, DoFType flag = All );

   void interpolateExtended( const std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                             const std::vector< VertexDoFFunction* >                                              srcFunctions,
                             uint_t                                                                               level,
                             DoFType                                                                              flag = All );

   /// assigns unique values to all data points
   /// this function is mainly used for petsc to get global identifier for all DoFs
   /// \tparam ValueType
   /// \param level
   void enumerate( uint_t level );

   // TODO: write more general version(s)
   ValueType getMaxValue( uint_t level, DoFType flag = All );
   ValueType getMinValue( uint_t level, DoFType flag = All );
   ValueType getMaxMagnitude( uint_t level, DoFType flag = All, bool mpiReduce = true );

   inline BoundaryCondition getBoundaryCondition() const { return boundaryCondition_; }

   template < typename SenderType, typename ReceiverType >
   inline void startCommunication( const uint_t& level ) const
   {
      if( isDummy() )
      {
         return;
      }
      communicators_.at( level )->template startCommunication< SenderType, ReceiverType >();
   }

   template < typename SenderType, typename ReceiverType >
   inline void endCommunication( const uint_t& level ) const
   {
      if( isDummy() )
      {
         return;
      }
      communicators_.at( level )->template endCommunication< SenderType, ReceiverType >();
   }

   template < typename SenderType, typename ReceiverType >
   inline void communicate( const uint_t& level ) const
   {
      if( isDummy() )
      {
         return;
      }
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

   void setLocalCommunicationMode( const communication::BufferedCommunicator::LocalCommunicationMode& localCommunicationMode );

   using Function< VertexDoFFunction< ValueType > >::isDummy;

 private:
   template < typename PrimitiveType >
   void interpolateByPrimitiveType( const ValueType& constant, uint_t level, DoFType flag = All ) const
   {
      if( isDummy() )
      {
         return;
      }
      this->startTiming( "Interpolate" );

      if( std::is_same< PrimitiveType, Vertex >::value )
      {
         for( const auto& it : this->getStorage()->getVertices() )
         {
            Vertex& vertex = *it.second;

            if( testFlag( boundaryCondition_.getBoundaryType( vertex.getMeshBoundaryFlag() ), flag ) )
            {
               vertexdof::macrovertex::interpolate( level, vertex, vertexDataID_, constant );
            }
         }
      } else if( std::is_same< PrimitiveType, Edge >::value )
      {
         for( const auto& it : this->getStorage()->getEdges() )
         {
            Edge& edge = *it.second;

            if( testFlag( boundaryCondition_.getBoundaryType( edge.getMeshBoundaryFlag() ), flag ) )
            {
               vertexdof::macroedge::interpolate( level, edge, edgeDataID_, constant );
            }
         }
      } else if( std::is_same< PrimitiveType, Face >::value )
      {
         for( const auto& it : this->getStorage()->getFaces() )
         {
            Face& face = *it.second;

            if( testFlag( boundaryCondition_.getBoundaryType( face.getMeshBoundaryFlag() ), flag ) )
            {
               vertexdof::macroface::interpolate( level, face, faceDataID_, constant );
            }
         }
      } else if( std::is_same< PrimitiveType, Cell >::value )
      {
         for( const auto& it : this->getStorage()->getCells() )
         {
            Cell& cell = *it.second;

            if( testFlag( boundaryCondition_.getBoundaryType( cell.getMeshBoundaryFlag() ), flag ) )
            {
               vertexdof::macrocell::interpolate( level, cell, cellDataID_, constant );
            }
         }
      }

      this->stopTiming( "Interpolate" );
   }

   inline void enumerate( uint_t level, ValueType& offset );

   using Function< VertexDoFFunction< ValueType > >::communicators_;
   using Function< VertexDoFFunction< ValueType > >::additiveCommunicators_;

   PrimitiveDataID< FunctionMemory< ValueType >, Vertex > vertexDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Edge >   edgeDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Face >   faceDataID_;
   PrimitiveDataID< FunctionMemory< ValueType >, Cell >   cellDataID_;

   BoundaryCondition boundaryCondition_;

   DoFType boundaryTypeToSkipDuringAdditiveCommunication_;

   /// friend P2Function for usage of enumerate
   friend class P2Function< ValueType >;
};

inline void projectMean( VertexDoFFunction< real_t >& pressure, VertexDoFFunction< real_t >& tmp, uint_t level )
{
   if( pressure.isDummy() )
   {
      return;
   }
   std::function< real_t( const hhg::Point3D& ) > ones = []( const hhg::Point3D& ) { return 1.0; };

   tmp.interpolate( ones, level );

   real_t numGlobalVertices = tmp.dotGlobal( tmp, level, hhg::All );
   real_t mean              = pressure.dotGlobal( tmp, level, hhg::All );

   pressure.assign( {1.0, -mean / numGlobalVertices}, {&pressure, &tmp}, level, hhg::All );
}

} // namespace vertexdof
} // namespace hhg
