cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)

find_package(openfpm 3.2.0 REQUIRED)
find_package(Threads)

if (openfpm_FOUND)
	message("OpenFPM found")
	

	add_executable(sph_dlb  main.cpp)
	get_property(OPENFPM_INCLUDES TARGET openfpm::binary_config PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
	get_property(OPENFPM_DEFINITION TARGET openfpm::binary_config PROPERTY INTERFACE_COMPILE_DEFINITIONS)
	get_property(OPENFPM_LIBS TARGET openfpm::binary_config PROPERTY INTERFACE_LINK_LIBRARIES)
	get_property(OPENFPM_COMPILE_OPTIONS TARGET openfpm::binary_config PROPERTY INTERFACE_COMPILE_OPTIONS)
	target_include_directories(sph_dlb PUBLIC ${OPENFPM_INCLUDES})
	target_compile_definitions(sph_dlb PUBLIC ${OPENFPM_DEFINITION})
	target_link_libraries(sph_dlb PUBLIC ${OPENFPM_LIBS})
	target_compile_options(sph_dlb PUBLIC ${OPENFPM_COMPILE_OPTIONS})
endif()

