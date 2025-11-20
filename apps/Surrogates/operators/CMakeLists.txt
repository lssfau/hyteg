add_library( opgen-mogli-diffusion
    P1ElementwiseDiffusion.cpp
    P1ElementwiseDiffusion.hpp
)
set_source_files_properties(
    P1ElementwiseDiffusion.cpp
    P1ElementwiseDiffusion.hpp
    PROPERTIES COMPILE_OPTIONS ${HYTEG_COMPILER_NATIVE_FLAGS}
)

if(HYTEG_BUILD_WITH_AVX AND WALBERLA_DOUBLE_ACCURACY)
    target_sources(opgen-mogli-diffusion PRIVATE
        P1ElementwiseDiffusion_cubes_const_vect_fused_quadloops_tab.cpp
        P1ElementwiseDiffusion_cubes_const_vect_fused_quadloops_tab.hpp
        P1ElementwiseDiffusion_IcosahedralShellMap_cubes_const_vect_fused_quadloops_tab.cpp
        P1ElementwiseDiffusion_IcosahedralShellMap_cubes_const_vect_fused_quadloops_tab.hpp
        P1ElementwiseDivKGrad_IcosahedralShellMap_cubes_const_vect_fused_quadloops_tab.cpp
        P1ElementwiseDivKGrad_IcosahedralShellMap_cubes_const_vect_fused_quadloops_tab.hpp
        P1ElementwiseDivKGrad_cubes_const_vect_fused_quadloops_tab.cpp
        P1ElementwiseDivKGrad_cubes_const_vect_fused_quadloops_tab.hpp
    )

    set_source_files_properties(
        P1ElementwiseDiffusion_cubes_const_vect_fused_quadloops_tab.cpp
        P1ElementwiseDiffusion_cubes_const_vect_fused_quadloops_tab.hpp
        P1ElementwiseDiffusion_IcosahedralShellMap_cubes_const_vect_fused_quadloops_tab.cpp
        P1ElementwiseDiffusion_IcosahedralShellMap_cubes_const_vect_fused_quadloops_tab.hpp
        P1ElementwiseDivKGrad_IcosahedralShellMap_cubes_const_vect_fused_quadloops_tab.cpp
        P1ElementwiseDivKGrad_IcosahedralShellMap_cubes_const_vect_fused_quadloops_tab.hpp
        P1ElementwiseDivKGrad_cubes_const_vect_fused_quadloops_tab.cpp
        P1ElementwiseDivKGrad_cubes_const_vect_fused_quadloops_tab.hpp
        PROPERTIES COMPILE_OPTIONS ${HYTEG_COMPILER_NATIVE_FLAGS}
    )
endif()

target_link_libraries( opgen-mogli-diffusion PRIVATE hyteg )