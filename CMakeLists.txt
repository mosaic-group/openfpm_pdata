cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)

if (POLICY CMP0074)
	cmake_policy(SET CMP0074 NEW)
endif()

set(openfpm_VERSION 5.0.0)

if (METIS_FOUND)
	set(DEFINE_HAVE_METIS "#define HAVE_METIS 1")
else()
	message( FATAL_ERROR "Metis is required in order to install OpenFPM" )
endif()

if (PARMETIS_FOUND)
	set(DEFINE_HAVE_PARMETIS "#define HAVE_PARMETIS 1")
else()
	message( FATAL_ERROR "ParMetis is required in order to install OpenFPM")
endif()

add_subdirectory (src)

file(READ ${CMAKE_CURRENT_SOURCE_DIR}/src/cmake/openfpmConfig-configure.cmake CMAKE_OPENFPM_CONFIG_VARS)
set(CMAKE_OPENFPM_CONFIG_VARS "${CMAKE_OPENFPM_CONFIG_VARS}\nmessage(STATUS \"Found OpenFPM version ${openfpm_VERSION} (\$\{CMAKE_CURRENT_LIST_FILE\})\")")
set(CMAKE_OPENFPM_CONFIG_VARS "${CMAKE_OPENFPM_CONFIG_VARS}\nset(OPENFPM_CUDA_ON_BACKEND \"${CUDA_ON_BACKEND}\")")
file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/src/cmake/openfpmConfig-Vars.cmake "${CMAKE_OPENFPM_CONFIG_VARS}")

get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
	set(DEFINE_HAVE_METIS ${DEFINE_HAVE_METIS} CACHE INTERNAL "")
	set(DEFINE_HAVE_PARMETIS ${DEFINE_HAVE_PARMETIS} CACHE INTERNAL "")
endif()

