
#pragma once

#include <core/all.h>

namespace hhg
{

template< typename SourceFunction, typename DestinationFunction >
class Operator
{
public:
  Operator(const std::shared_ptr<PrimitiveStorage> & storage, uint_t minLevel, uint_t maxLevel)
    : storage_(storage), minLevel_(minLevel), maxLevel_(maxLevel)
  {
  }

  virtual ~Operator()
  {
  }

  void apply( SourceFunction& src, DestinationFunction& dst, size_t level, DoFType flag, UpdateType updateType = Replace );

  void smooth_gs( DestinationFunction& dst, SourceFunction& rhs, size_t level, DoFType flag );

 protected:

  virtual void apply_impl( SourceFunction& src, DestinationFunction& dst, size_t level, DoFType flag, UpdateType updateType = Replace ) = 0;
  virtual void smooth_gs_impl( DestinationFunction& dst, SourceFunction& rhs, size_t level, DoFType flag ) = 0;

  const std::shared_ptr< PrimitiveStorage > storage_;
  const uint_t minLevel_;
  const uint_t maxLevel_;
};


template< typename SourceFunction, typename DestinationFunction >
void Operator< SourceFunction, DestinationFunction  >::apply( SourceFunction& src, DestinationFunction& dst, size_t level, DoFType flag, UpdateType updateType )
{
  apply_impl( src, dst, level, flag, updateType );
}

template< typename SourceFunction, typename DestinationFunction >
void Operator< SourceFunction, DestinationFunction  >::smooth_gs( DestinationFunction& dst, SourceFunction& rhs, size_t level, DoFType flag )
{
  smooth_gs_impl( dst, rhs, level, flag );
}

}
