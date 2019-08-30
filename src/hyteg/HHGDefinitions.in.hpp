#pragma once

#cmakedefine HHG_BUILD_WITH_PETSC
#cmakedefine HHG_BUILD_WITH_EIGEN
#cmakedefine HHG_USE_GENERATED_KERNELS
#cmakedefine HHG_P1_COLORING

#ifdef HHG_USE_GENERATED_KERNELS
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

#ifdef HHG_P1_COLORING
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