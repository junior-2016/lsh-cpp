cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(lsh_cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# find boost library
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost REQUIRED)

# find GSL
find_package(GSL REQUIRED)

# require: sudo apt-get install python-matplotlib python-numpy python-tk python2.7-dev
find_package(PythonLibs 2.7 REQUIRED)

# ArchLinux OS don't link numpy headers to /usr/include/python2.7.
set(NUMPY_INCLUDE "/usr/lib/python2.7/site-packages/numpy/core/include/")

# find eigen & intel_mkl
find_package(Eigen3 3.3.7 REQUIRED)
if ("$ENV{MKLROOT}" STREQUAL "")
    message(FATAL_ERROR "You don't hava intel mkl.")
endif ()
message("MKLROOT = $ENV{MKLROOT}")
set(MKLROOT "$ENV{MKLROOT}")
add_definitions(-DCPU_THREAD_NUMBERS=8) # specify cpu thread numbers for mkl

# find stmpct include directories
set(LIBSTMPCT_PUBLIC_HDR "third-party/streaming-percentiles/cpp/src/include/stmpct/")

# source file path
file(GLOB Self_include_file "src/include/*.h")
file(GLOB Test_file "src/test/*.cpp" "src/test/*.h")
file(GLOB Benchmark_file "src/benchmark/*.h" "src/benchmark/*.cpp")

# option for test and benchmark
OPTION(TEST "build test" ON)
OPTION(BENCHMARK "build benchmark" ON)
set(TEST ON)
set(BENCHMARK ON)

# 有关GCC优化相关的更多选项见:
# https://stackoverflow.com/questions/14492436/g-optimization-beyond-o3-ofast
# https://gcc.gnu.org/onlinedocs/gcc-9.1.0/gcc/Optimize-Options.html
# for release
set(CMAKE_CXX_FLAGS "-Ofast -DNDEBUG")
# for debug
set(CMAKE_CXX_FLAGS_DEBUG "-O0")

set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Wextra -Wconversion -pedantic -Werror=return-type")

# 有关GCC对本地指令集启用和优化的选项见: https://gcc.gnu.org/onlinedocs/gcc-9.1.0/gcc/x86-Options.html
# (-march=native会开启本机所有可行的指令集,比如sse/avx等; -mtune=native会根据这些指令集范围进行优化)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} -fopenmp -DMKL_LP64 -m64 -mavx2 -msse4.2 -mfma")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${GCC_COVERAGE_COMPILE_FLAGS} -fopenmp -DMKL_LP64 -m64 -mavx2 -msse4.2 -mfma")

# use libc++ library if using clang++ compiler.
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_definitions(-DUSE_CLANG_COMPILER) # 添加 USE_CLANG_COMPILER 宏定义.注意添加宏的语法是 -Dmacro.
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++") # 使用 llvm 的 libc++ 标准库
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -stdlib=libc++")
endif ()

# check gcc version
OPTION(USE_CXX_PARALLISM_TS "use c++17 parallism TS" OFF)
if (CMAKE_COMPILER_IS_GNUCXX AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 9.1)
    # only gcc 9.1.0 can use c++17 parallism ts
    add_definitions("-DUSE_CXX_PARALLISM_TS")
    set(USE_CXX_PARALLISM_TS ON)
    # find intel tbb as gcc parallism TS backend.
    if ("$ENV{TBBROOT}" STREQUAL "")
        message(FATAL_ERROR "You don't have intel tbb.")
    endif ()
    set(TBBROOT "$ENV{TBBROOT}")
    message("TBBROOT = ${TBBROOT}")
    find_library(
            TBBLIB
            NAMES libtbb.so.2 libtbbmalloc.so.2 libtbbmalloc_proxy.so.2 libtbb_preview.so.2
            PATHS "${TBBROOT}/lib/intel64/gcc4.7"
    )
    message(STATUS "TBBLIB = " ${TBBLIB})
else ()
    set(USE_CXX_PARALLISM_TS OFF)
endif ()
# unset cmake option variable cache,see:https://stackoverflow.com/questions/22481647/cmake-if-else-with-option
unset(USE_CXX_PARALLISM_TS CACHE)

# executable setting
include_directories(${Boost_INCLUDE_DIRS} ${NUMPY_INCLUDE} ${MKLROOT}/include ${EIGEN3_INCLUDE_DIR}
        ${LIBSTMPCT_PUBLIC_HDR} ${TBBROOT}/include)

# libraries setting
set(LSH_CPP_REQUIRE_LIBS ${Boost_LIBRARIES}
        ${PYTHON_LIBRARIES}
        GSL::gsl GSL::gslcblas
        ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a
        -Wl,--start-group
        ${MKLROOT}/lib/intel64/libmkl_cdft_core.a
        ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a
        ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a
        ${MKLROOT}/lib/intel64/libmkl_core.a
        ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a
        -Wl,--end-group
        -lgomp -lpthread -lm -ldl)
if (USE_CXX_PARALLISM_TS)
    set(LSH_CPP_REQUIRE_LIBS ${LSH_CPP_REQUIRE_LIBS} ${TBBLIB})
endif ()

if (TEST)
    add_executable(lsh_cpp_test ${Self_include_file} ${Test_file})
    target_include_directories(lsh_cpp_test PRIVATE ${PYTHON_INCLUDE_DIRS})
    target_link_libraries(lsh_cpp_test ${LSH_CPP_REQUIRE_LIBS})
ENDIF ()
if (BENCHMARK)
    add_executable(lsh_cpp_benchmark ${Self_include_file} ${Benchmark_file})
    target_include_directories(lsh_cpp_benchmark PRIVATE ${PYTHON_INCLUDE_DIRS})
    target_link_libraries(lsh_cpp_benchmark ${LSH_CPP_REQUIRE_LIBS})
ENDIF ()

