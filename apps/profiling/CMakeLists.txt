    add_executable       ( ProfilingApp ProfilingApp.cpp )
    target_link_libraries( ProfilingApp hyteg walberla::core constant_stencil_operator )