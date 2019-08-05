#pragma once

#include "tinyhhg_core/forms/Form.hpp"

namespace hhg {

class P1Form : public Form
{
 public:
   virtual ~P1Form() {}

   // 2D P1
   virtual void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const { WALBERLA_ABORT( "Not implemented." ); }

   // 3D P1
   virtual void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const { WALBERLA_ABORT( "Not implemented." ); }
};

} // namespace hhg
