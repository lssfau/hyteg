message(VERBOSE "Detecting AVX support...")

if(CMAKE_CXX_COMPILER_ID MATCHES "Intel" OR CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
  set(HYTEG_COMPILER_NATIVE_FLAGS "-xhost")
else()
  set(HYTEG_COMPILER_NATIVE_FLAGS "-march=native")
endif()

try_compile(HYTEG_PLATFORM_SUPPORTS_AVX
   SOURCE_FROM_CONTENT main.cpp "
      #include <immintrin.h>
      int main() {
         const __m256d result = _mm256_mul_pd(_mm256_set_pd(4.0,5.0,6.0,7.0), _mm256_set_pd(5.0,6.0,7.0,8.0));
      }
      "
   COMPILE_DEFINITIONS ${HYTEG_COMPILER_NATIVE_FLAGS}
   OUTPUT_VARIABLE OUT
)

message(DEBUG ${OUT})

if(${HYTEG_PLATFORM_SUPPORTS_AVX})
  message(STATUS "Platform supports AVX.")
else()
  message(STATUS "Platform does NOT support AVX.")
endif()
