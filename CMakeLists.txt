cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata VERSION 5.2.0 LANGUAGES C CXX)

if (POLICY CMP0074)
	cmake_policy(SET CMP0074 NEW)
endif()

if (NOT openfpm_VERSION)
	set(openfpm_VERSION "${PROJECT_VERSION}")
endif()

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

set(OPENFPM_PACKAGE_CONFIG_BUILD_DIR
	"${CMAKE_CURRENT_BINARY_DIR}/openfpm_cmake")
set(OPENFPM_CONFIG_VARS_FILE
	"${OPENFPM_PACKAGE_CONFIG_BUILD_DIR}/openfpmConfig-Vars.cmake")
set(OPENFPM_CONFIG_FILE
	"${OPENFPM_PACKAGE_CONFIG_BUILD_DIR}/openfpmConfig.cmake")
file(MAKE_DIRECTORY "${OPENFPM_PACKAGE_CONFIG_BUILD_DIR}")
configure_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/src/cmake/openfpm-config.cmake.in"
	"${OPENFPM_CONFIG_FILE}"
	COPYONLY)

add_subdirectory (src)

unset(CMAKE_OPENFPM_CONFIG_VARS CACHE)
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/src/cmake/openfpmConfig-configure.cmake
	OPENFPM_CONFIG_VARS_CONTENT)
set(OPENFPM_CONFIG_VARS_CONTENT
	"if(NOT OPENFPM_CONFIG_INITIALIZED)\n${OPENFPM_CONFIG_VARS_CONTENT}")
set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\nset(OPENFPM_CUDA_ON_BACKEND \"${CUDA_ON_BACKEND}\")")
set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\ninclude(CMakeFindDependencyMacro)")
set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\nif(NOT Boost_DIR)\n  set(Boost_DIR \"${Boost_DIR}\")\nendif()")
set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\nfind_dependency(Boost 1.72 COMPONENTS context fiber filesystem iostreams program_options)")
set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\nif(NOT MPI_CXX_COMPILER)\n  set(MPI_CXX_COMPILER \"${MPI_CXX_COMPILER}\")\nendif()")
set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\nfind_dependency(MPI COMPONENTS CXX)")
set(_openfpm_runtime_library_dirs
	"${CMAKE_INSTALL_PREFIX}/openfpm_devices/lib"
	"${CMAKE_INSTALL_PREFIX}/openfpm_vcluster/lib"
	"${CMAKE_INSTALL_PREFIX}/openfpm_pdata/lib")
set(_openfpm_runtime_libraries
	${MPI_C_LIBRARIES}
	${MPI_CXX_LIBRARIES}
	${PARMETIS_LIBRARIES}
	${HDF5_LIBRARIES}
	${PETSC_LIBRARIES}
	${SuiteSparse_LIBRARIES}
	${Boost_LIBRARIES})
set(_openfpm_mpi_root "")
if (IS_ABSOLUTE "${MPI_CXX_COMPILER}")
	get_filename_component(_openfpm_mpi_bindir "${MPI_CXX_COMPILER}" DIRECTORY)
	get_filename_component(_openfpm_mpi_root "${_openfpm_mpi_bindir}" DIRECTORY)
endif()
foreach(_openfpm_runtime_library IN LISTS _openfpm_runtime_libraries)
	if (IS_ABSOLUTE "${_openfpm_runtime_library}" AND
		EXISTS "${_openfpm_runtime_library}" AND
		NOT IS_DIRECTORY "${_openfpm_runtime_library}")
		get_filename_component(_openfpm_runtime_library_dir
			"${_openfpm_runtime_library}" DIRECTORY)
		list(APPEND _openfpm_runtime_library_dirs
			"${_openfpm_runtime_library_dir}")
	endif()
endforeach()
foreach(_openfpm_runtime_root IN ITEMS
	"${_openfpm_mpi_root}" "${MPI_HOME}" "${METIS_ROOT}" "${PARMETIS_ROOT}" "${HDF5_ROOT}"
	"${PETSC_ROOT}" "${SUITESPARSE_ROOT}" "${OPENBLAS_ROOT}"
	"${LIBHILBERT_ROOT}" "${BOOST_ROOT}")
	if (_openfpm_runtime_root)
		foreach(_openfpm_runtime_libdir IN ITEMS lib lib64)
			if (IS_DIRECTORY
				"${_openfpm_runtime_root}/${_openfpm_runtime_libdir}")
				list(APPEND _openfpm_runtime_library_dirs
					"${_openfpm_runtime_root}/${_openfpm_runtime_libdir}")
			endif()
		endforeach()
	endif()
endforeach()
list(APPEND _openfpm_runtime_library_dirs
	${Boost_LIBRARY_DIRS} ${LIBHILBERT_LIBRARY_DIRS})
list(REMOVE_DUPLICATES _openfpm_runtime_library_dirs)
set(OPENFPM_CONFIG_VARS_CONTENT
	"${OPENFPM_CONFIG_VARS_CONTENT}\nset(OPENFPM_MPIEXEC_EXECUTABLE \"${MPIEXEC_EXECUTABLE}\")")
set(OPENFPM_CONFIG_VARS_CONTENT
	"${OPENFPM_CONFIG_VARS_CONTENT}\nset(OPENFPM_RUNTIME_LIBRARY_DIRS \"${_openfpm_runtime_library_dirs}\")")
set(OPENFPM_CONFIG_VARS_CONTENT
	"${OPENFPM_CONFIG_VARS_CONTENT}\ninclude(\"\${CMAKE_CURRENT_LIST_DIR}/OpenFPMRun.cmake\")")
if (CUDA_ON_BACKEND STREQUAL "OpenMP")
	set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\nfind_dependency(OpenMP COMPONENTS CXX)")
endif()
if (CUDA_ON_BACKEND STREQUAL "METAL")
	set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\nset(OPENFPM_MOLTENVK_HIP_CLANG \"${OPENFPM_MOLTENVK_HIP_CLANG}\")")
	set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\nset(OPENFPM_MOLTENVK_LLVM_OPT \"${OPENFPM_MOLTENVK_LLVM_OPT}\")")
	set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\nset(OPENFPM_MOLTENVK_LLVM_LINK \"${OPENFPM_MOLTENVK_LLVM_LINK}\")")
	set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\nset(OPENFPM_MOLTENVK_CLSPV \"${OPENFPM_MOLTENVK_CLSPV}\")")
	set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\nset(OPENFPM_MOLTENVK_SPIRV_VAL \"${OPENFPM_MOLTENVK_SPIRV_VAL}\")")
	set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\nset(OPENFPM_CHIPSTAR_ROOT \"${OPENFPM_CHIPSTAR_ROOT}\")")
	set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\nset(OPENFPM_CHIPSTAR_GENERATED_INCLUDE \"${OPENFPM_CHIPSTAR_GENERATED_INCLUDE}\")")
	set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\ninclude(\"\${CMAKE_CURRENT_LIST_DIR}/../openfpm_devices/cmake/OpenFPMMoltenVK.cmake\")")
endif()
set(OPENFPM_CONFIG_VARS_CONTENT "${OPENFPM_CONFIG_VARS_CONTENT}\ninclude(\"\${CMAKE_CURRENT_LIST_DIR}/../openfpm_devices/cmake/OpenFPMGPU.cmake\")")
set(OPENFPM_CONFIG_VARS_CONTENT
	"${OPENFPM_CONFIG_VARS_CONTENT}\nset(OPENFPM_CONFIG_INITIALIZED TRUE)\nendif()")

set(OPENFPM_NUMERICS_BUILT FALSE)
if (ENABLE_NUMERICS)
	set(OPENFPM_NUMERICS_BUILT TRUE)
endif()
set(OPENFPM_PETSC_FOUND_AT_BUILD FALSE)
if (PETSC_FOUND)
	set(OPENFPM_PETSC_FOUND_AT_BUILD TRUE)
endif()
configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/src/cmake/openfpm-config-components.cmake.in
	${OPENFPM_PACKAGE_CONFIG_BUILD_DIR}/openfpmConfig-components.cmake
	@ONLY)
file(READ
	${OPENFPM_PACKAGE_CONFIG_BUILD_DIR}/openfpmConfig-components.cmake
	OPENFPM_COMPONENT_CONFIG)
set(OPENFPM_CONFIG_VARS_CONTENT
	"${OPENFPM_CONFIG_VARS_CONTENT}\n${OPENFPM_COMPONENT_CONFIG}")
file(WRITE "${OPENFPM_CONFIG_VARS_FILE}" "${OPENFPM_CONFIG_VARS_CONTENT}")

get_directory_property(hasParent PARENT_DIRECTORY)
if(hasParent)
	set(DEFINE_HAVE_METIS ${DEFINE_HAVE_METIS} CACHE INTERNAL "")
	set(DEFINE_HAVE_PARMETIS ${DEFINE_HAVE_PARMETIS} CACHE INTERNAL "")
endif()
