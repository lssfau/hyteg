#pragma once

#cmakedefine HYTEG_BUILD_WITH_PETSC
#cmakedefine HYTEG_BUILD_WITH_EIGEN
#cmakedefine HYTEG_USE_GENERATED_KERNELS
#cmakedefine HYTEG_P1_COLORING

#ifdef HYTEG_USE_GENERATED_KERNELS
namespace hyteg {
namespace globalDefines {
constexpr bool useGeneratedKernels = true;
} // namespace globalDefines
} // namespace hyteg
#else
namespace hyteg {
namespace globalDefines {
constexpr bool useGeneratedKernels = false;
} // namesapce globalDefines
} // namespace hyteg
#endif

#ifdef HYTEG_P1_COLORING
namespace hyteg {
namespace globalDefines {
constexpr bool useP1Coloring = true;
} // namespace globalDefines
} // namespace hyteg
#else
namespace hyteg {
namespace globalDefines {
constexpr bool useP1Coloring = false;
} // namesapce globalDefines
} // namespace hyteg
#endif