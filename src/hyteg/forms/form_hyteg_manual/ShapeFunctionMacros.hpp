
// =====================================
//  P1 shape functions on unit triangle
// =====================================
#ifdef DEFINE_P1_SHAPE_FUNCTIONS_TRIANGLE
#define SF_M0 L1
#define SF_M1 L2
#define SF_M2 L3
#endif

#ifdef UNDEFINE_P1_SHAPE_FUNCTIONS_TRIANGLE
#undef SF_M0
#undef SF_M1
#undef SF_M2
#endif

// ========================================
//  P1 shape functions on unit tetrahedron
// ========================================
#ifdef DEFINE_P1_SHAPE_FUNCTIONS_TET
#define SF_M0 L1
#define SF_M1 L2
#define SF_M2 L3
#define SF_M3 L4
#endif

#ifdef UNDEFINE_P1_SHAPE_FUNCTIONS_TET
#undef SF_M0
#undef SF_M1
#undef SF_M2
#undef SF_M3
#endif

// =====================================
//  P2 shape functions on unit triangle
// =====================================
#ifdef DEFINE_P2_SHAPE_FUNCTIONS_TRIANGLE
#define SF_N0 ( L1 * ( 2.0 * L1 - 1.0 ) )
#define SF_N1 ( L2 * ( 2.0 * L2 - 1.0 ) )
#define SF_N2 ( L3 * ( 2.0 * L3 - 1.0 ) )
#define SF_N3 ( 4.0 * L2 * L3 )
#define SF_N4 ( 4.0 * L1 * L3 )
#define SF_N5 ( 4.0 * L1 * L2 )
#endif

#ifdef UNDEFINE_P2_SHAPE_FUNCTIONS_TRIANGLE
#undef SF_N0
#undef SF_N1
#undef SF_N2
#undef SF_N3
#undef SF_N4
#undef SF_N5
#endif

// ========================================
//  P2 shape functions on unit tetrahedron
// ========================================
#ifdef DEFINE_P2_SHAPE_FUNCTIONS_TET
#define SF_N0 ( L1 * ( 2.0 * L1 - 1.0 ) )
#define SF_N1 ( L2 * ( 2.0 * L2 - 1.0 ) )
#define SF_N2 ( L3 * ( 2.0 * L3 - 1.0 ) )
#define SF_N3 ( L4 * ( 2.0 * L4 - 1.0 ) )
#define SF_N4 ( 4.0 * L3 * L4 )
#define SF_N5 ( 4.0 * L2 * L4 )
#define SF_N6 ( 4.0 * L2 * L3 )
#define SF_N7 ( 4.0 * L1 * L4 )
#define SF_N8 ( 4.0 * L1 * L3 )
#define SF_N9 ( 4.0 * L1 * L2 )
#endif

#ifdef UNDEFINE_P2_SHAPE_FUNCTIONS_TET
#undef SF_N0
#undef SF_N1
#undef SF_N2
#undef SF_N3
#undef SF_N4
#undef SF_N5
#undef SF_N6
#undef SF_N7
#undef SF_N8
#undef SF_N9
#endif

// ====================================================
//  Derivatives of P2 shape functions on unit triangle
// ====================================================

#ifdef DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE

#define DxiN0 ( 1.0 - 4.0 * L1 )
#define DetaN0 ( 1.0 - 4.0 * L1 )

#define DxiN1 ( 4.0 * L2 - 1.0 )
#define DetaN1 ( 0.0 )

#define DxiN2 ( 0.0 )
#define DetaN2 ( 4.0 * L3 - 1.0 )

#define DxiN3 ( 4.0 * L3 )
#define DetaN3 ( 4.0 * L2 )

#define DxiN4 ( -4.0 * L3 )
#define DetaN4 ( 4.0 * ( L1 - L3 ) )

#define DxiN5 ( 4.0 * ( L1 - L2 ) )
#define DetaN5 ( -4.0 * L2 )

#endif

#ifdef UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TRIANGLE
#undef DetaN0
#undef DetaN1
#undef DetaN2
#undef DetaN3
#undef DetaN4
#undef DetaN5
#undef DxiN0
#undef DxiN1
#undef DxiN2
#undef DxiN3
#undef DxiN4
#undef DxiN5
#endif

// =======================================================
//  Derivatives of P2 shape functions on unit tetrahedron
// =======================================================
#ifdef DEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET

#define DX1N0 ( L2 + L3 + L4 - 3.0 * L1 )
#define DX2N0 ( L2 + L3 + L4 - 3.0 * L1 )
#define DX3N0 ( L2 + L3 + L4 - 3.0 * L1 )

#define DX1N1 ( 4.0 * L2 - 1.0 )
#define DX2N1 ( 0.0 )
#define DX3N1 ( 0.0 )

#define DX1N2 ( 0.0 )
#define DX2N2 ( 4.0 * L3 - 1.0 )
#define DX3N2 ( 0.0 )

#define DX1N3 ( 0.0 )
#define DX2N3 ( 0.0 )
#define DX3N3 ( 4.0 * L4 - 1.0 )

#define DX1N4 ( 0.0 )
#define DX2N4 ( 4.0 * L4 )
#define DX3N4 ( 4.0 * L3 )

#define DX1N5 ( 4.0 * L4 )
#define DX2N5 ( 0.0 )
#define DX3N5 ( 4.0 * L2 )

#define DX1N6 ( 4.0 * L3 )
#define DX2N6 ( 4.0 * L2 )
#define DX3N6 ( 0.0 )

#define DX1N7 ( - 4.0 * L4 )
#define DX2N7 ( - 4.0 * L4 )
#define DX3N7 ( 4.0 * ( L1 - L4 ) )

#define DX1N8 ( - 4.0 * L3 )
#define DX2N8 ( 4.0 * ( L1 - L3 ) )
#define DX3N8 ( - 4.0 * L3 )

#define DX1N9 ( 4.0 * ( L1 - L2 ) )
#define DX2N9 ( - 4.0 * L2 )
#define DX3N9 ( - 4.0 * L2 )

#endif

#ifdef UNDEFINE_P2_SHAPE_FUNCTION_DERIVATIVES_TET

#undef DX1N0
#undef DX2N0
#undef DX3N0

#undef DX1N1
#undef DX2N1
#undef DX3N1

#undef DX1N2
#undef DX2N2
#undef DX3N2

#undef DX1N3
#undef DX2N3
#undef DX3N3

#undef DX1N4
#undef DX2N4
#undef DX3N4

#undef DX1N5
#undef DX2N5
#undef DX3N5

#undef DX1N6
#undef DX2N6
#undef DX3N6

#undef DX1N7
#undef DX2N7
#undef DX3N7

#undef DX1N8
#undef DX2N8
#undef DX3N8

#undef DX1N9
#undef DX2N9
#undef DX3N9

#endif
