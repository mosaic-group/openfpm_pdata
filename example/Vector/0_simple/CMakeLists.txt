cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)

find_package(openfpm 3.2.0 REQUIRED)
find_package(Threads)

if (openfpm_FOUND)
	message("OpenFPM found")
	

	add_executable(vect  main.cpp)
	get_property(OPENFPM_INCLUDES TARGET openfpm::binary_config PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
	get_property(OPENFPM_DEFINITION TARGET openfpm::binary_config PROPERTY INTERFACE_COMPILE_DEFINITIONS)
	get_property(OPENFPM_LIBS TARGET openfpm::binary_config PROPERTY INTERFACE_LINK_LIBRARIES)
	target_include_directories(vect PUBLIC ${OPENFPM_INCLUDES})
	target_compile_definitions(vect PUBLIC ${OPENFPM_DEFINITION})
	target_link_libraries(vect PUBLIC ${OPENFPM_LIBS})
endif()

