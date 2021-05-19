cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)

if (POLICY CMP0074)
	cmake_policy(SET CMP0074 OLD)
endif ()

enable_testing()

if (POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

set(openfpm_VERSION 3.3.0)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake_modules/)

set(BOOST_INCLUDE ${Boost_INCLUDE_DIR} CACHE PATH "Include directory for BOOST")
set(PETSC_ROOT CACHE PATH "If compiling with linear algebra indicate the PETSC root directory")
set(PARMETIS_ROOT CACHE PATH "Parmetis root directory")
set(METIS_ROOT CACHE PATH "Metis root directory")
set(LIBHILBERT_ROOT CACHE PATH "LibHilbert root path")
set(HDF5_ROOT CACHE PATH "HDF5 root path")
set(EIGEN3_ROOT CACHE PATH "Eigen3 include path")
set(LIBHILBERT_ROOT CACHE PATH "LibHilbert root path")
set(SUITESPARSE_ROOT CACHE PATH "The suitesparse root path")
set(TINYOBJLOADER_ROOT CACHE PATH "TinyObjLoader library path")
set(SE_CLASS1 CACHE BOOL "Activate compilation with SE_CLASS1")
set(SE_CLASS2 CACHE BOOL "Activate compilation with SE_CLASS2")
set(SE_CLASS3 CACHE BOOL "Activate compilation with SE_CLASS3")
set(ACTION_ON_ERROR CACHE STRING "Action to perform in case of error")
set(PROFILE_WITH_SCOREP CACHE BOOL "Enable profiling with scorep")
set(ENV{PETSC_DIR} ${PETSC_ROOT})
#set(ENV{HDF5_ROOT} ${HDF5_ROOT})
set(ENV{EIGEN3_ROOT} ${EIGEN3_ROOT})
set(METIS_DIR ${METIS_ROOT})
set(PARMETIS_DIR ${PARMETIS_ROOT})
set(OPENBLAS_ROOT CACHE PATH "Root path for blas library")
set(SuiteSparse_ROOT ${SUITESPARSE_ROOT})
set(CUDA_ON_CPU CACHE BOOL "Make Cuda work on heap")
set(CPACK_RUN_INSTALL_DEPENDENCIES CACHE BOOL "Set to true if we are creating deb or RPM packages")
set(ENABLE_GARBAGE_INJECTOR CACHE BOOL "Enable the injector of garbage in the memory allocator")
set(ENABLE_VCLUSTER_GARBAGE_INJECTOR CACHE BOOL "Enable the injector of garbage in the vcluster memory buffers")
set(HIP_ENABLE CACHE BOOL "Enable HIP compiler")
set(AMD_ARCH_COMPILE "gfx900" CACHE STRING "AMD gpu architecture used to compile kernels")

# Enabling real GPU is stronger than using CUDA_ON_CPU
if (ENABLE_GPU)
	set(CUDA_ON_CPU OFF)
	# Test CLang
	if (NOT HIP_ENABLE)
		find_package(CUDA)
		set(CMAKE_CUDA_COMPILER_LIBRARY_ROOT ${CUDA_TOOLKIT_ROOT_DIR})
        	enable_language(CUDA)
	endif()
endif()

set (CMAKE_CXX_STANDARD 14)
set (CMAKE_CUDA_STANDARD 14)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3")

if (OPENBLAS_ROOT)
	set(ENV{LD_LIBRARY_PATH} "$ENV{LD_LIBRARY_PATH}:${OPENBLAS_ROOT}/lib")
	set(ENV{DYLD_LIBRARY_PATH} "$ENV{DYLD_LIBRARY_PATH}:${OPENBLAS_ROOT}/lib")
endif()


set(ENV{PATH} "$ENV{PATH}:${HDF5_ROOT}/bin")
set(HDF5_PREFER_PARALLEL TRUE)

set(Vc_DIR "${Vc_ROOT}/lib/cmake/Vc/")
message("Searching Vc in ${Vc_DIR}")

find_package(Boost 1.72.0 COMPONENTS unit_test_framework iostreams program_options system filesystem OPTIONAL_COMPONENTS fiber context)
find_package(MPI)
find_package(PETSc)
find_package(HDF5)
find_package(Eigen3)
find_package(LibHilbert)
find_package(METIS)
find_package(ParMETIS)
find_package(TinyObjLoader )
find_package(BLAS)
find_package(LAPACK)
find_package(Eigen3)
find_package(SuiteSparse OPTIONAL_COMPONENTS UMFPACK)
find_package(Vc)
find_package(OpenMP)
find_package(HIP)

set(CMAKE_SKIP_BUILD_RPATH TRUE)

if(HIP_FOUND)
        set(DEFINE_HIP_GPU "#define HIP_GPU")
        set(DEFINE_CUDIFY_USE_HIP "#define CUDIFY_USE_HIP")
	file(WRITE hip_enabled 1)
else()
	file(WRITE hip_enabled 0)
endif()

if(HIP_FOUND)
        set(DEFINE_CUDA_GPU "#define CUDA_GPU")
endif()

if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

if(PROFILE_WITH_SCOREP)
	set(CMAKE_CXX_COMPILER_LAUNCHER "scorep")
	set(CMAKE_CC_COMPILER_LAUNCHER "scorep")
	set(CMAKE_CUDA_COMPILER_LAUNCHER "scorep")
endif()


if(ENABLE_GPU AND (CUDA_FOUND OR HIP_FOUND))
	set(OPENFPM_INIT_FILE "initialize/initialize_wrapper_cuda.cu")
else()
	set(OPENFPM_INIT_FILE "initialize/initialize_wrapper_cpu.cpp")
endif()

###### CONFIG.h FILE ######

if(SE_CLASS1)
	set(DEFINE_SE_CLASS1 "#define SE_CLASS1")
endif()

if(ACTION_ON_ERROR)
	set(DEFINE_ACTION_ON_ERROR "#define ${ACTION_ON_ERROR}")
endif()

if(SE_CLASS2)
	set(DEFINE_SE_CLASS2 "#define SE_CLASS2")
endif()

if(SE_CLASS3)
	set(DEFINE_SE_CLASS3 "#define SE_CLASS3")
endif()

if(PETSC_FOUND)
	set(DEFINE_HAVE_PETSC "#define HAVE_PETSC")
endif()

if(SCAN_COVERTY)
	set(DEFINE_SCAN_COVERTY "#define COVERTY_SCAN")
endif()

if(TEST_PERFORMANCE)
	set(DEFINE_PERFORMANCE_TEST "#define PERFORMANCE_TEST")
endif()

if (METIS_FOUND)
	set(DEFINE_HAVE_METIS "#define HAVE_METIS 1")
else()
	file(WRITE error_code "201")
	message( FATAL_ERROR "Metis is required in order to install OpenFPM" )
endif()

if (PARMETIS_FOUND)
	set(DEFINE_HAVE_PARMETIS "#define HAVE_PARMETIS 1")
else()
	file(WRITE error_code "208")
	message( FATAL_ERROR "ParMetis is required in order to install OpenFPM")
endif()


if(MPI_FOUND)
	get_filename_component(OPENFPM_MPI_DEP "${MPI_C_INCLUDE_DIRS}" DIRECTORY)
	file(READ ${CMAKE_SOURCE_DIR}/src/cmake/openfpmConfig-configure.cmake CMAKE_OPENFPM_CONFIG_VARS)
	file(WRITE ${CMAKE_SOURCE_DIR}/src/cmake/openfpmConfigVars-configure.cmake "${CMAKE_OPENFPM_CONFIG_VARS}\nset(CMAKE_PREFIX_PATH ${OPENFPM_MPI_DEP}/)")
	set(DEFINE_HAVE_MPI "#define HAVE_MPI")
else()
	file(WRITE error_code "200")
	message( FATAL_ERROR "MPI is required in order to install OpenFPM" )
endif()

if (Boost_FOUND)
        set(DEFINE_HAVE_BOOST "#define HAVE_BOOST")
        set(DEFINE_HAVE_BOOST_IOSTREAMS "#define HAVE_BOOST_IOSTREAMS")
        set(DEFINE_HAVE_BOOST_PROGRAM_OPTIONS "#define HAVE_BOOST_PROGRAM_OPTIONS")
        set(DEFINE_HAVE_BOOST_UNIT_TEST_FRAMEWORK "#define HAVE_BOOST_UNIT_TEST_FRAMEWORK")
        if (Boost_CONTEXT_FOUND)
                set(DEFINE_HAVE_BOOST_CONTEXT "#define HAVE_BOOST_CONTEXT")
		set(OPTIONAL_BOOST_LIBS "-lboost_context")
        else()
                #if context is not there CUDA_ON_CPU cannot be activated
                set(CUDA_ON_CPU OFF)
        endif()
        if (Boost_FIBER_FOUND)
                set(DEFINE_HAVE_BOOST_FIBER "#define HAVE_BOOST_FIBER")
		string(CONCAT OPTIONAL_BOOST_LIBS ${OPTIONAL_BOOST_LIBS}  " -lboost_fiber")
        endif()
	file(WRITE optional_boost_libs "${OPTIONAL_BOOST_LIBS}")
else()
	file(WRITE error_code "202")
	message( FATAL_ERROR "BOOST is required in order to install OpenFPM" )
endif()

if(ENABLE_GPU AND CUDA_FOUND)
        set(DEFINE_CUDA_GPU "#define CUDA_GPU")
endif()

if(CUDA_ON_CPU)
        set(DEFINE_CUDA_GPU "#define CUDA_GPU")
endif()

if(HDF5_FOUND)
	if (HDF5_IS_PARALLEL)
		set(DEFINE_HAVE_HDF5 "#define HAVE_HDF5")
	else()
		file(WRITE error_code "207")
		message( FATAL_ERROR "HDF5 found ${HDF5_INCLUDE_DIRS} does not have parallel support, OpenFPM require it" )
	endif()
else()
	file(WRITE error_code "207")
	message( FATAL_ERROR "HDF5 with parallel support is required in order to install OpenFPM" )
endif()

if(EIGEN3_FOUND)
    message("-- EIGEN3 found.")
    set(DEFINE_HAVE_EIGEN "#define HAVE_EIGEN")
else()
    message("-- EIGEN3 not found!")
endif()

if(LIBHILBERT_FOUND)
	set(DEFINE_HAVE_LIBHILBERT "#define HAVE_LIBHILBERT 1")
else()
	file(WRITE error_code "210")
	message( FATAL_ERROR "LibHilbert is required in order to install OpenFPM")
endif()

if(SUITESPARSE_FOUND AND SuiteSparse_UMFPACK_FOUND)
	set(DEFINE_HAVE_SUITESPARSE "#define HAVE_SUITESPARSE")
endif()

if (Vc_FOUND)
	set(DEFINE_HAVE_VCDEVEL)
else()
	file(WRITE error_code "211")
	message( FATAL_ERROR "Vc is required in roder to install OpenFPM")
endif()

if (ENABLE_GARBAGE_INJECTOR)
	set(DEFINE_GARBAGE_INJECTOR "#define GARBAGE_INJECTOR")
endif()

if (ENABLE_VCLUSTER_GARBAGE_INJECTOR)
        set(DEFINE_VCLUSTER_GARBAGE_INJECTOR "#define VCLUSTER_GARBAGE_INJECTOR")
endif()

if(APPLE)
	set(DEFINE_HAVE_OSX "#define HAVE_OSX")
endif()

if(TINYOBJLOADER_FOUND)
	set(DEFINE_HAVE_TINYOBJLOADER "#define HAVE_TINYOBJLOADER 1")
endif()

file(WRITE error_code "0")
file(WRITE cuda_lib "${CUDA_cudart_static_LIBRARY} ${CUDA_cudadevrt_LIBRARY}")
if(CUDA_ON_CPU)
	file(WRITE cuda_include "-I${CUDA_INCLUDE_DIRS} -D__NVCC__ -DCUDART_VERSION=11000 -DCUDA_ON_CPU")
	file(WRITE cuda_on_cpu "YES")
else()
	file(WRITE cuda_include "-I${CUDA_INCLUDE_DIRS}")
	file(WRITE cuda_on_cpu "NO")
endif()
file(WRITE mpi_include "-I${MPI_C_INCLUDE_DIRS}")
file(WRITE mpi_libs "${MPI_C_LINK_FLAGS} ${MPI_C_LIBRARIES}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/config/config_cmake.h.in ${CMAKE_CURRENT_SOURCE_DIR}/src/config/config.h)

add_subdirectory (openfpm_io)
add_subdirectory (openfpm_numerics)

file(WRITE cuda_options "${WARNING_SUPPRESSION_AND_OPTION_NVCC_TEXT}")

add_subdirectory (src)

#################### CPack to create auto installing packages

include(InstallRequiredSystemLibraries)

string(REPLACE "." ";" VERSION_LIST ${openfpm_VERSION})
list(GET VERSION_LIST 0 OPENFPM_VERSION_MAJOR)
list(GET VERSION_LIST 1 OPENFPM_VERSION_MINOR)
list(GET VERSION_LIST 2 OPENFPM_VERSION_PATCH)

if (CPACK_RUN_INSTALL_DEPENDENCIES)

	set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "OpenFPM distributed data-structures")
	set(CPACK_PACKAGE_VENDOR "IBirdSoft")
	set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
	set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/gpl-3.0.txt")
	set(CPACK_PACKAGE_VERSION_MAJOR ${OPENFPM_VERSION_MAJOR})
	set(CPACK_PACKAGE_VERSION_MINOR ${OPENFPM_VERSION_MINOR})
	set(CPACK_PACKAGE_VERSION_PATCH ${OPENFPM_VERSION_PATCH})
	set(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")
	set(CPACK_PACKAGE_INSTALL_DIRECTORY /usr/local/openfpm)
	set(CPACK_PACKAGING_INSTALL_PREFIX /usr/local/openfpm)
	set(CPACK_PACKAGE_HOMEPAGE_URL http://openfpm.mpi-cbg.de)
	set(CPACK_RPM_PACKAGE_AUTOREQPROV NO)
	set(CPACK_PACKAGE_CONTACT incardon@mpi-cbg.de)
	set(CPACK_DEBIAN_PACKAGE_MAINTAINER Pietro Incardona)
	set(CPACK_RPM_POST_INSTALL_SCRIPT_FILE ${CMAKE_CURRENT_SOURCE_DIR}/src/scripts/postinst)
	set(CPACK_POSTFLIGHT_OPENFPM_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/src/scripts/postflight)
	set(CPACK_DEBIAN_PACKAGE_CONTROL_EXTRA "${CMAKE_CURRENT_SOURCE_DIR}/src/script/postinst")
	set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")

	install(FILES $ENV{DEP_PACKING}/openfpm_vars
		DESTINATION source
              	COMPONENT OpenFPM)

        install(DIRECTORY $ENV{DEP_PACKING}/BOOST
                DESTINATION dependencies/
                COMPONENT OpenFPM)

        install(DIRECTORY $ENV{DEP_PACKING}/EIGEN
                DESTINATION dependencies/
                COMPONENT OpenFPM)

        install(DIRECTORY $ENV{DEP_PACKING}/HDF5
                DESTINATION dependencies/
                COMPONENT OpenFPM)

        install(DIRECTORY $ENV{DEP_PACKING}/LIBHILBERT
                DESTINATION dependencies/
                COMPONENT OpenFPM)

        install(DIRECTORY $ENV{DEP_PACKING}/METIS
                DESTINATION dependencies/
                COMPONENT OpenFPM)

        install(DIRECTORY $ENV{DEP_PACKING}/MPI
                DESTINATION dependencies/
                COMPONENT OpenFPM)

        install(DIRECTORY $ENV{DEP_PACKING}/OPENBLAS
                DESTINATION dependencies/
                COMPONENT OpenFPM)

        install(DIRECTORY $ENV{DEP_PACKING}/PARMETIS
                DESTINATION dependencies/
                COMPONENT OpenFPM)

        install(DIRECTORY $ENV{DEP_PACKING}/PETSC
                DESTINATION dependencies/
                COMPONENT OpenFPM)

        install(DIRECTORY $ENV{DEP_PACKING}/SUITESPARSE
                DESTINATION dependencies/
                COMPONENT OpenFPM)

	install(DIRECTORY $ENV{DEP_PACKING}/VCDEVEL
                DESTINATION dependencies/
                COMPONENT OpenFPM)

endif()

include(CPack)

cpack_add_component(OpenFPM
                    DISPLAY_NAME OpenFPM
                    DESCRITION OpenFPM binary files)


