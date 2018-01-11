#pragma once

#include "core/DataTypes.h"

namespace hhg {

using walberla::real_c;

template<typename ValueType>
class P2Function : public Function< P2Function< ValueType > >
{
public:

  P2Function( const std::string& name, const std::shared_ptr< PrimitiveStorage > & storage, uint_t minLevel, uint_t maxLevel ) :
      Function< P2Function< ValueType > >( name, storage, minLevel, maxLevel ),
      vertexDoFFunction_( std::make_shared< vertexdof::VertexDoFFunction< ValueType > >( name + "_VertexDoF", storage, minLevel, maxLevel ) ),
      edgeDoFFunction_(   std::make_shared<              EdgeDoFFunction< ValueType > >( name + "_VertexDoF", storage, minLevel, maxLevel ) )
  {}

  std::shared_ptr< vertexdof::VertexDoFFunction< ValueType > > getVertexDoFFunction() const { return vertexDoFFunction_; }
  std::shared_ptr<              EdgeDoFFunction< ValueType > > getEdgeDoFFunction()   const { return edgeDoFFunction_;   }

private:

    inline void
    interpolate_impl( std::function< ValueType( const Point3D&, const std::vector< ValueType >& ) >& expr,
                      const std::vector< P2Function< ValueType > * > srcFunctions, uint_t level, DoFType flag = All )
    {
      std::vector< vertexdof::VertexDoFFunction< real_t > * > vertexDoFFunctions;
      std::vector<              EdgeDoFFunction< real_t > * > edgeDoFFunctions;

      for ( const auto & function : srcFunctions )
      {
        vertexDoFFunctions.push_back( function->vertexDoFFunction_.get() );
        edgeDoFFunctions.push_back  ( function->edgeDoFFunction_.get() );
      }

      vertexDoFFunction_->interpolateExtended( expr, vertexDoFFunctions, level, flag );
      edgeDoFFunction_->interpolateExtended( expr, edgeDoFFunctions, level, flag );
    }

    inline void
    assign_impl( const std::vector< ValueType > scalars, const std::vector< P2Function< ValueType >* > functions,
                 uint_t level, DoFType flag = All )
    {
      std::vector< vertexdof::VertexDoFFunction< real_t > * > vertexDoFFunctions;
      std::vector<              EdgeDoFFunction< real_t > * > edgeDoFFunctions;

      for ( const auto & function : functions )
      {
        vertexDoFFunctions.push_back( function->vertexDoFFunction_.get() );
        edgeDoFFunctions.push_back  ( function->edgeDoFFunction_.get() );
      }

      vertexDoFFunction_->assign( scalars, vertexDoFFunctions, level, flag );
      edgeDoFFunction_->assign  ( scalars, edgeDoFFunctions,   level, flag );
    }

    inline void
    add_impl( const std::vector< ValueType > scalars, const std::vector< P2Function< ValueType >* > functions,
              uint_t level, DoFType flag = All )
    {
      std::vector< vertexdof::VertexDoFFunction< real_t > * > vertexDoFFunctions;
      std::vector<              EdgeDoFFunction< real_t > * > edgeDoFFunctions;

      for ( const auto & function : functions )
      {
        vertexDoFFunctions.push_back( function->vertexDoFFunction_.get() );
        edgeDoFFunctions.push_back  ( function->edgeDoFFunction_.get() );
      }

      vertexDoFFunction_->add( scalars, vertexDoFFunctions, level, flag );
      edgeDoFFunction_->add  ( scalars, edgeDoFFunctions,   level, flag );
    }

    inline real_t
    dot_impl( P2Function< ValueType >& rhs, uint_t level, DoFType flag = All )
    {
      real_t sum = real_c( 0 );
      sum += vertexDoFFunction_->dot( *rhs.vertexDoFFunction_, level, flag);
      sum += edgeDoFFunction_->dot  ( *rhs.edgeDoFFunction_,   level, flag);
      return sum;
    }

    inline void
    prolongate_impl( uint_t sourceLevel, DoFType flag = All )
    {
      WALBERLA_ABORT( "P2Function - Prolongate not implemented!" );
    }

    inline void
    prolongateQuadratic_impl( uint_t sourceLevel, DoFType flag = All )
    {
      WALBERLA_ABORT( "P2Function - Prolongate (quadratic) not implemented!" );
    }

    inline void
    restrict_impl( uint_t sourceLevel, DoFType flag = All )
    {
      WALBERLA_ABORT( "P2Function - Restrict not implemented!" );
    }

    inline void
    enumerate_impl( uint_t level, uint_t& num )
    {
      vertexDoFFunction_->enumerate( level, num );
      edgeDoFFunction_->enumerate  ( level, num );
    }

  std::shared_ptr< vertexdof::VertexDoFFunction< ValueType > > vertexDoFFunction_;
  std::shared_ptr< EdgeDoFFunction< ValueType > >              edgeDoFFunction_;

};

}
