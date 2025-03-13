add_library( opgen-mogli-diffusion
    P1ElementwiseDiffusion.cpp
    P1ElementwiseDiffusion.hpp
    P1ElementwiseDiffusion_cubes_const_vect_polycse.cpp
    P1ElementwiseDiffusion_cubes_const_vect_polycse.hpp
)
set_source_files_properties(
    P1ElementwiseDiffusion.cpp
    P1ElementwiseDiffusion.hpp
    P1ElementwiseDiffusion_cubes_const_vect_polycse.cpp
    P1ElementwiseDiffusion_cubes_const_vect_polycse.hpp
    PROPERTIES COMPILE_OPTIONS ${HYTEG_COMPILER_NATIVE_FLAGS}
)
target_link_libraries( opgen-mogli-diffusion PRIVATE hyteg )