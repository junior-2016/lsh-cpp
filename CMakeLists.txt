cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(lsh_cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

# 有关GCC优化相关的更多选项见:
# https://stackoverflow.com/questions/14492436/g-optimization-beyond-o3-ofast
# https://gcc.gnu.org/onlinedocs/gcc-9.1.0/gcc/Optimize-Options.html
# for release
set(CMAKE_CXX_FLAGS "-Ofast")
# for debug
set(CMAKE_CXX_FLAGS_DEBUG "-O0")

set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Wextra -Wconversion -pedantic -Werror=return-type")

# 有关GCC对本地指令集启用和优化的选项见: https://gcc.gnu.org/onlinedocs/gcc-9.1.0/gcc/x86-Options.html
# (-march=native会开启本机所有可行的指令集,比如sse/avx等; -mtune=native会根据这些指令集范围进行优化)
# 可以考虑用 -march=native -mtune=native -fopenmp 将全部可用指令集打开;
# 但如果sse/avx全开可能导致两种指令混用,最后性能严重下降,所以这里仅启用 sse4.2 (前提是本机允许sse4.2指令集)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} -msse4.2 -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${GCC_COVERAGE_COMPILE_FLAGS} -msse4.2 -fopenmp")

# find boost library
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED)

# find GSL
find_package(GSL REQUIRED)

# find anaconda python library and numpy headers
set(PYTHON_INCLUDE_PATH "$ENV{HOME}/anaconda3/include/python3.7m")
set(NUMPY_CONDA_INCLUDE "$ENV{HOME}/anaconda3/lib/python3.7/site-packages/numpy/core/include")
set(PYTHON_LIBRARIES "$ENV{HOME}/anaconda3/lib/libpython3.7m.so")

# add google abseil-cpp library dependency
add_subdirectory(./third-party/abseil-cpp)

# source file path
set(SELF_INCLUDE src/include/)
file(GLOB Test_file "src/test/*.cpp" "src/test/*.h")
file(GLOB Benchmark_file "src/benchmark/*.h" "src/benchmark/*.cpp")

OPTION(USE_SIMD "use simd or not" OFF)
OPTION(TEST "build test" ON)
OPTION(BENCHMARK "build benchmark" OFF)
set(TEST ON)
set(BENCHMARK ON)
set(USE_SIMD OFF)

IF (USE_SIMD)
    ADD_DEFINITIONS(-DUSE_SIMD)
    find_package(xsimd REQUIRED CONFIG) # find xsimd package
ENDIF (USE_SIMD)

if (Boost_FOUND AND GSL_FOUND)
    include_directories(${Boost_INCLUDE_DIRS} ${NUMPY_CONDA_INCLUDE} ${SELF_INCLUDE})
    if (TEST)
        add_executable(lsh_cpp_test ${Test_file})
        target_include_directories(lsh_cpp_test PRIVATE ${PYTHON_INCLUDE_PATH})
        target_link_libraries(lsh_cpp_test absl::hash ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}
                GSL::gsl GSL::gslcblas)
    ENDIF ()
    if (BENCHMARK)
        add_executable(lsh_cpp_benchmark ${Benchmark_file})
        target_include_directories(lsh_cpp_benchmark PRIVATE ${PYTHON_INCLUDE_PATH})
        target_link_libraries(lsh_cpp_benchmark absl::hash ${Boost_LIBRARIES} ${PYTHON_LIBRARIES}
                GSL::gsl GSL::gslcblas)
    ENDIF ()
endif ()
