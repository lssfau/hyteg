

if (HYTEG_BUILD_WITH_PETSC)
add_executable       ( MultiSinkerBenchmark MultiSinkerBenchmark.cpp )
target_link_libraries( MultiSinkerBenchmark hyteg walberla::core mixed_operator )

add_executable       ( SolViBenchmark SolViBenchmark.cpp )
target_link_libraries( SolViBenchmark hyteg walberla::core mixed_operator )

add_executable       ( EGSmoothViscTests EGSmoothViscTests.cpp )
target_link_libraries( EGSmoothViscTests hyteg walberla::core mixed_operator )

add_executable       ( EGJumpTests.cpp EGJumpTests.cpp )
target_link_libraries( EGJumpTests.cpp hyteg walberla::core mixed_operator )
endif ()
