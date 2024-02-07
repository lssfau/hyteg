message(VERBOSE "Detecting AVX support...")

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel" OR CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
  set(HYTEG_COMPILER_NATIVE_FLAGS "-xhost")
else()
  set(HYTEG_COMPILER_NATIVE_FLAGS "-march=native")
endif()

try_compile(HYTEG_PLATFORM_SUPPORTS_AVX
   "${CMAKE_CURRENT_BINARY_DIR}"
   "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/TestAvxSupport.cpp"
   COMPILE_DEFINITIONS ${HYTEG_COMPILER_NATIVE_FLAGS}
   OUTPUT_VARIABLE OUT
)

message(DEBUG ${OUT})

if(${HYTEG_PLATFORM_SUPPORTS_AVX})
  message(STATUS "Platform supports AVX.")
else()
  message(STATUS "Platform does NOT support AVX.")
endif()
