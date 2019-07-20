cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
project(lsh_cpp)

set(CMAKE_CXX_STANDARD 17)

set(GCC_COVERAGE_COMPILE_FLAGS "-Wall -Wextra -Wconversion -pedantic -Werror=return-type")

# for release :
set(CMAKE_CXX_FLAGS "-O3")
# for debug :
set(CMAKE_CXX_FLAGS_DEBUG "-O0")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${GCC_COVERAGE_COMPILE_FLAGS}")

# find boost library
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.65.1)

file(GLOB source_files "src/include/*.h" "src/include/*.cpp")

set(PYTHON_INCLUDE_PATH "$ENV{HOME}/anaconda3/include/python3.7m")
set(NUMPY_CONDA_INCLUDE "$ENV{HOME}/anaconda3/lib/python3.7/site-packages/numpy/core/include")
set(PYTHON_LIBRARIES "$ENV{HOME}/anaconda3/lib/libpython3.7m.so")

if (Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    include_directories(${NUMPY_CONDA_INCLUDE})
    add_executable(lsh_cpp src/main.cpp ${source_files})
    target_include_directories(lsh_cpp PRIVATE ${PYTHON_INCLUDE_PATH})
    target_link_libraries(lsh_cpp ${Boost_LIBRARIES})
    target_link_libraries(lsh_cpp ${PYTHON_LIBRARIES})
endif ()

