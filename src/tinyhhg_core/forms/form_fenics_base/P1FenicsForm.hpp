#pragma once

#include "tinyhhg_core/forms/P1Form.hpp"

#include "tinyhhg_core/fenics/fenics.hpp"
#include "tinyhhg_core/fenics/ufc_traits.hpp"

#ifdef _MSC_VER
#  pragma warning(push, 0)
#endif

// P1
#include "tinyhhg_core/forms/form_fenics_generated/p1_diffusion.h"
#include "tinyhhg_core/forms/form_fenics_generated/p1_div.h"
#include "tinyhhg_core/forms/form_fenics_generated/p1_div_K_grad.h"
#include "tinyhhg_core/forms/form_fenics_generated/p1_divt.h"
#include "tinyhhg_core/forms/form_fenics_generated/p1_mass.h"
#include "tinyhhg_core/forms/form_fenics_generated/p1_polar_laplacian.h"
#include "tinyhhg_core/forms/form_fenics_generated/p1_polar_mass.h"
#include "tinyhhg_core/forms/form_fenics_generated/p1_pspg.h"
#include "tinyhhg_core/forms/form_fenics_generated/p1_stokes_epsilon.h"
#include "tinyhhg_core/forms/form_fenics_generated/p1_tet_diffusion.h"
#include "tinyhhg_core/forms/form_fenics_generated/p1_tet_div_tet.h"
#include "tinyhhg_core/forms/form_fenics_generated/p1_tet_divt_tet.h"
#include "tinyhhg_core/forms/form_fenics_generated/p1_tet_mass.h"
#include "tinyhhg_core/forms/form_fenics_generated/p1_tet_pspg_tet.h"


// P2
#include "tinyhhg_core/forms/form_fenics_generated/p2_diffusion.h"
#include "tinyhhg_core/forms/form_fenics_generated/p2_div.h"
#include "tinyhhg_core/forms/form_fenics_generated/p2_divt.h"
#include "tinyhhg_core/forms/form_fenics_generated/p2_mass.h"
#include "tinyhhg_core/forms/form_fenics_generated/p2_tet_diffusion.h"
#include "tinyhhg_core/forms/form_fenics_generated/p2_tet_div_tet.h"
#include "tinyhhg_core/forms/form_fenics_generated/p2_tet_divt_tet.h"
#include "tinyhhg_core/forms/form_fenics_generated/p2_tet_mass.h"
#include "tinyhhg_core/forms/form_fenics_generated/p2_tet_pspg_tet.h"

// P1 to P2
#include "tinyhhg_core/forms/form_fenics_generated/p1_to_p2_divt.h"
#include "tinyhhg_core/forms/form_fenics_generated/p1_to_p2_tet_divt_tet.h"

// P2 to P1
#include "tinyhhg_core/forms/form_fenics_generated/p2_to_p1_div.h"
#include "tinyhhg_core/forms/form_fenics_generated/p2_to_p1_tet_div_tet.h"

#ifdef _MSC_VER
#  pragma warning(pop)
#endif

namespace hyteg {

template < class UFCOperator2D, class UFCOperator3D = fenics::UndefinedAssembly >
class P1FenicsForm : public P1Form
{
 public:
   void integrate( const std::array< Point3D, 3 >& coords, Point3D& out ) const override
   {
      Matrix3r localStiffnessMatrix;
      real_t   fenicsCoords[6];
      fenicsCoords[0] = coords[0][0];
      fenicsCoords[1] = coords[0][1];
      fenicsCoords[2] = coords[1][0];
      fenicsCoords[3] = coords[1][1];
      fenicsCoords[4] = coords[2][0];
      fenicsCoords[5] = coords[2][1];
      UFCOperator2D gen;
      gen.tabulate_tensor( localStiffnessMatrix.data(), nullptr, fenicsCoords, 0 );
      out[0] = localStiffnessMatrix( 0, 0 );
      out[1] = localStiffnessMatrix( 0, 1 );
      out[2] = localStiffnessMatrix( 0, 2 );
   }

   void integrate( const std::array< Point3D, 4 >& coords, Point4D& out ) const override
   {
      typename fenics::UFCTrait< UFCOperator3D >::LocalStiffnessMatrix_T localStiffnessMatrix;

      // Flattening the offset array to be able to pass it to the fenics routines.
      double geometricOffsetsArray[12];
      for ( uint_t cellVertex = 0; cellVertex < 4; cellVertex++ )
      {
         for ( int coordinate = 0; coordinate < 3; coordinate++ )
         {
            geometricOffsetsArray[cellVertex * 3 + uint_c( coordinate )] = coords[cellVertex][coordinate];
         }
      }

      UFCOperator3D gen;
      gen.tabulate_tensor( localStiffnessMatrix.data(), NULL, geometricOffsetsArray, 0 );

      out[0] = localStiffnessMatrix( 0, 0 );
      out[1] = localStiffnessMatrix( 0, 1 );
      out[2] = localStiffnessMatrix( 0, 2 );
      out[3] = localStiffnessMatrix( 0, 3 );
   }

   bool assemble2D() const override { return !std::is_same< UFCOperator2D, hyteg::fenics::NoAssemble >::value; }

   bool assemble3D() const override { return !std::is_same< UFCOperator3D, hyteg::fenics::NoAssemble >::value; }

   bool assembly2DDefined() const override { return !std::is_same< UFCOperator2D, hyteg::fenics::UndefinedAssembly >::value; }

   bool assembly3DDefined() const override { return !std::is_same< UFCOperator3D, hyteg::fenics::UndefinedAssembly >::value; }
};

} // namespace hyteg
