#pragma once

#include "tinyhhg_core/geometry/GeometryMap.hpp"

namespace hyteg {

class Form
{
 public:
   virtual ~Form() {}

   virtual bool assemble2D() const = 0;

   virtual bool assemble3D() const = 0;

   virtual bool assembly2DDefined() const = 0;

   virtual bool assembly3DDefined() const = 0;

   std::shared_ptr< GeometryMap > geometryMap;
};

} // namespace hyteg