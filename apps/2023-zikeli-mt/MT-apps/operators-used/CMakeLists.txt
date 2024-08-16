add_compile_options( "-Wno-shadow" )

if (WALBERLA_BUILD_WITH_HALF_PRECISION_SUPPORT)
    add_library( opgen-mt-diffusion

        P1ElementwiseDiffusion_cubes_const_float16.cpp
        P1ElementwiseDiffusion_cubes_const_float16.hpp
        P1ElementwiseDiffusion_cubes_const_float32.cpp
        P1ElementwiseDiffusion_cubes_const_float32.hpp
        P1ElementwiseDiffusion_cubes_const_float64.cpp
        P1ElementwiseDiffusion_cubes_const_float64.hpp

        P1ElementwiseDiffusion_cubes_const_vect_float32.cpp
        P1ElementwiseDiffusion_cubes_const_vect_float32.hpp
        P1ElementwiseDiffusion_cubes_const_vect_float64.cpp
        P1ElementwiseDiffusion_cubes_const_vect_float64.hpp
    )

    set_source_files_properties(

        P1ElementwiseDiffusion_cubes_const_vect_float32.cpp
        P1ElementwiseDiffusion_cubes_const_vect_float32.hpp
        P1ElementwiseDiffusion_cubes_const_vect_float64.cpp
        P1ElementwiseDiffusion_cubes_const_vect_float64.hpp

        PROPERTIES COMPILE_OPTIONS ${HYTEG_COMPILER_NATIVE_FLAGS}
    )

    target_compile_features(opgen-mt-diffusion PUBLIC cxx_std_23) # probably not necessary since this is inherited by module core.
else ()
    add_library( opgen-mt-diffusion

        P1ElementwiseDiffusion_cubes_const_float32.cpp
        P1ElementwiseDiffusion_cubes_const_float32.hpp
        P1ElementwiseDiffusion_cubes_const_float64.cpp
        P1ElementwiseDiffusion_cubes_const_float64.hpp

        P1ElementwiseDiffusion_cubes_const_vect_float32.cpp
        P1ElementwiseDiffusion_cubes_const_vect_float32.hpp
        P1ElementwiseDiffusion_cubes_const_vect_float64.cpp
        P1ElementwiseDiffusion_cubes_const_vect_float64.hpp
    )

    set_source_files_properties(

        P1ElementwiseDiffusion_cubes_const_vect_float32.cpp
        P1ElementwiseDiffusion_cubes_const_vect_float32.hpp
        P1ElementwiseDiffusion_cubes_const_vect_float64.cpp
        P1ElementwiseDiffusion_cubes_const_vect_float64.hpp

        PROPERTIES COMPILE_OPTIONS ${HYTEG_COMPILER_NATIVE_FLAGS}
    )
endif ()
