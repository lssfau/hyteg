#pragma once

#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/fenics/ufc_traits.hpp"
#include "tinyhhg_core/p2functionspace/P2Form.hpp"

// P1
#include "tinyhhg_core/p1functionspace/generated/p1_diffusion.h"
#include "tinyhhg_core/p1functionspace/generated/p1_div.h"
#include "tinyhhg_core/p1functionspace/generated/p1_div_K_grad.h"
#include "tinyhhg_core/p1functionspace/generated/p1_divt.h"
#include "tinyhhg_core/p1functionspace/generated/p1_mass.h"
#include "tinyhhg_core/p1functionspace/generated/p1_polar_laplacian.h"
#include "tinyhhg_core/p1functionspace/generated/p1_polar_mass.h"
#include "tinyhhg_core/p1functionspace/generated/p1_pspg.h"
#include "tinyhhg_core/p1functionspace/generated/p1_stokes_epsilon.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_diffusion.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_div_tet.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_divt_tet.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_mass.h"
#include "tinyhhg_core/p1functionspace/generated/p1_tet_pspg_tet.h"

// P2
#include "tinyhhg_core/p2functionspace/generated/p2_diffusion.h"
#include "tinyhhg_core/p2functionspace/generated/p2_div.h"
#include "tinyhhg_core/p2functionspace/generated/p2_divt.h"
#include "tinyhhg_core/p2functionspace/generated/p2_mass.h"
#include "tinyhhg_core/p2functionspace/generated/p2_pspg.h"
#include "tinyhhg_core/p2functionspace/generated/p2_tet_diffusion.h"
#include "tinyhhg_core/p2functionspace/generated/p2_tet_div_tet.h"
#include "tinyhhg_core/p2functionspace/generated/p2_tet_divt_tet.h"
#include "tinyhhg_core/p2functionspace/generated/p2_tet_mass.h"
#include "tinyhhg_core/p2functionspace/generated/p2_tet_pspg_tet.h"

// P1 to P2
#include "tinyhhg_core/mixedoperators/generated/p1_to_p2_divt.h"
#include "tinyhhg_core/mixedoperators/generated/p1_to_p2_tet_divt_tet.h"

// P2 to P1
#include "tinyhhg_core/mixedoperators/generated/p2_to_p1_div.h"
#include "tinyhhg_core/mixedoperators/generated/p2_to_p1_tet_div_tet.h"

namespace hhg {

template < class UFCOperator2D, class UFCOperator3D = fenics::UndefinedAssembly >
class P2FenicsForm : public P2Form
{
 public:

   // ---------------------------
   //  2D versions for triangles
   // ---------------------------
   void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const override
   {
      Matrix6r localStiffnessMatrix;
      computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 0, 0 );
      out[1] = localStiffnessMatrix( 0, 1 );
      out[2] = localStiffnessMatrix( 0, 2 );
   }

   void integrateEdgeToVertex( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      Matrix6r localStiffnessMatrix;
      computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 0, 3 );
      out[1] = localStiffnessMatrix( 0, 4 );
      out[2] = localStiffnessMatrix( 0, 5 );
   }

   void integrateVertexToEdge( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      Matrix6r localStiffnessMatrix;
      computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 5, 0 );
      out[1] = localStiffnessMatrix( 5, 1 );
      out[2] = localStiffnessMatrix( 5, 2 );
   }

   void integrateEdgeToEdge( const std::array< Point3D, 3 >& coords, Point3D& out ) const
   {
      Matrix6r localStiffnessMatrix;
      computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 5, 3 );
      out[1] = localStiffnessMatrix( 5, 4 );
      out[2] = localStiffnessMatrix( 5, 5 );
   }

   // ----------------------------
   //  3D versions for tetrahedra
   // ----------------------------
   void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const override {
      Matrix10r localStiffnessMatrix;
      computeLocalStiffnessMatrix( coords, localStiffnessMatrix );
      out[0] = localStiffnessMatrix( 0, 0 );
      out[1] = localStiffnessMatrix( 0, 1 );
      out[2] = localStiffnessMatrix( 0, 2 );
      out[3] = localStiffnessMatrix( 0, 3 );
   }

   void integrateEdgeToVertex( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2FenicsForm::integrateEdgeToVertex() not implemented for 3D!" );
   }

   void integrateVertexToEdge( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2FenicsForm::integrateVertexToEdge() not implemented for 3D!" );
   }

   void integrateEdgeToEdge( const std::array< Point3D, 4 >& coords, Point4D& out ) const
   {
      WALBERLA_ABORT( "P2FenicsForm::integrateEdgeToEdge() not implemented for 3D!" );
   }

   // -------------

   bool assemble2D() const override { return !std::is_same< UFCOperator2D, hhg::fenics::NoAssemble >::value; }

   bool assemble3D() const override { return !std::is_same< UFCOperator3D, hhg::fenics::NoAssemble >::value; }

   bool assembly2DDefined() const override { return !std::is_same< UFCOperator2D, hhg::fenics::UndefinedAssembly >::value; }

   bool assembly3DDefined() const override { return !std::is_same< UFCOperator3D, hhg::fenics::UndefinedAssembly >::value; }

 private:
   void computeLocalStiffnessMatrix( const std::array< Point3D, 3 >& coords, Matrix6r& localStiffnessMatrix ) const
   {
      real_t fenicsCoords[6];
      fenicsCoords[0] = coords[0][0];
      fenicsCoords[1] = coords[0][1];
      fenicsCoords[2] = coords[1][0];
      fenicsCoords[3] = coords[1][1];
      fenicsCoords[4] = coords[2][0];
      fenicsCoords[5] = coords[2][1];
      UFCOperator2D gen;
      gen.tabulate_tensor( localStiffnessMatrix.data(), nullptr, fenicsCoords, 0 );
   }

   void computeLocalStiffnessMatrix( const std::array< Point3D, 4 >& coords, Matrix10r& localStiffnessMatrix ) const
   {
      real_t fenicsCoords[12];
      for( int node = 0; node < 4; ++node ) {
        for( int dim = 0; dim < 3; ++dim ) {
          fenicsCoords[node*3+dim] = coords[node][dim];
        }
      }
      UFCOperator3D gen;
      gen.tabulate_tensor( localStiffnessMatrix.data(), nullptr, fenicsCoords, 0 );
   }
};

} // namespace hhg
