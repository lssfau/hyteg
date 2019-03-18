#pragma once

#include "tinyhhg_core/geometry/GeometryMap.hpp"

namespace hhg {

class Form
{
 public:
   virtual ~Form() {}

   // 2D P1
   virtual void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const { WALBERLA_ABORT( "Not implemented." ); }

   // 3D P1
   virtual void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const { WALBERLA_ABORT( "Not implemented." ); }

   virtual bool assemble2D() const = 0;

   virtual bool assemble3D() const = 0;

   virtual bool assembly2DDefined() const = 0;

   virtual bool assembly3DDefined() const = 0;

   std::shared_ptr< GeometryMap > geometryMap;
};

} // namespace hhg