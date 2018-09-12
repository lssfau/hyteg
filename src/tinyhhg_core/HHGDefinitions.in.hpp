#pragma once

#cmakedefine HHG_BUILD_WITH_PETSC
#cmakedefine HHG_BUILD_WITH_EIGEN
#cmakedefine HHG_USE_GENERATED_KERNELS

#ifdef HHG_USE_GENERATED_KERNELS
namespace hhg {
namespace globalDefines {
constexpr bool useGeneratedKernels = true;
} // namespace globalDefines
} // namespace hhg
#else
namespace hhg {
namespace globalDefines {
constexpr bool useGeneratedKernels = false;
} // namesapce globalDefines
} // namespace hhg
#endif