add_executable(lbm d2q9.c main.c)
target_include_directories(lbm PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(lbm PRIVATE m)
set_property(TARGET lbm
             PROPERTY C_STANDARD 11)

find_package(OpenMP)
if(OpenMP_C_FOUND)
  target_link_libraries(lbm PRIVATE OpenMP::OpenMP_C)
endif()

# Compile Options
include(compile-flags-helpers)
include(${PROJECT_SOURCE_DIR}/optimization_flags.cmake)

if (DEFINED ADDITIONAL_BENCHMARK_COMPILE_OPTIONS)
  add_compiler_option_to_target_type(lbm Benchmark PRIVATE ${ADDITIONAL_BENCHMARK_COMPILE_OPTIONS})
endif()

foreach(compile_type IN ITEMS Release RelWithDebInfo)
  add_compiler_option_to_target_type(lbm ${compile_type} PRIVATE ${ADDITIONAL_RELEASE_COMPILE_OPTIONS})
  add_linker_option_to_target_type(lbm ${compile_type} PRIVATE ${ADDITIONAL_RELEASE_LINK_OPTIONS})
endforeach()

add_compiler_option_to_target_type(lbm Debug PRIVATE ${ADDITIONAL_DEBUG_COMPILE_OPTIONS})

# Linker Options

if (DEFINED ADDITIONAL_BENCHMARK_LINK_OPTIONS)
  add_linker_option_to_target_type(lbm Benchmark PRIVATE ${ADDITIONAL_BENCHMARK_LINK_OPTIONS})
endif()

add_sanitizers_to_target(lbm Debug PRIVATE address undefined)

include(CheckIPOSupported)
check_ipo_supported(RESULT result)
if((result) AND USE_IPO)
  set_property(TARGET lbm PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

install(TARGETS lbm RUNTIME DESTINATION bin)