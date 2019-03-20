#pragma once

#include "tinyhhg_core/form/Form.hpp"

namespace hhg {

class P2Form : public Form
{
 public:
   virtual ~P2Form() {}

   // 2D P2 VertexDoF
   virtual void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const { WALBERLA_ABORT( "Not implemented." ); }

   // 3D P2 VertexDoF
   virtual void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const { WALBERLA_ABORT( "Not implemented." ); }
};

} // namespace hhg