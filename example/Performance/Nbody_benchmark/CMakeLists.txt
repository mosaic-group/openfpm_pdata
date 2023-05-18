# Kokkos minimally requires 3.10 right now,
# but your project can set it higher
cmake_minimum_required(VERSION 3.10)

# Projects can safely mix languages - must have C++ support
# Kokkos flags will only apply to C++ files
project(Example CXX)

# You need this for using Kokkos_ROOT variable
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.12.0")
  message(STATUS "Setting policy CMP0074 to use <Package>_ROOT variables")
  cmake_policy(SET CMP0074 NEW)
endif()

find_package(OpenMP)
find_package(Threads)
find_package(MPI)
find_package(Boost 1.75.0 COMPONENTS program_options iostreams filesystem fiber context REQUIRED)

# Look for an installed Kokkos
find_package(openfpm 4.1.0 REQUIRED)

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3")

add_definitions( -D__NVCC__ -DCUDART_VERSION=11000 )
add_compile_options(  -ftree-vectorize  -fstrict-aliasing -ffast-math  -march=skylake -mavx2 -mtune=skylake -mprefer-vector-width=256  )

# -fopt-info-vec-missed

add_executable(nbody main.cpp)

# This is the only thing required to set up compiler/linker flags
target_link_libraries(nbody  openfpm::binary_config)
target_link_libraries(nbody OpenMP::OpenMP_CXX)
