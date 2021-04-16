cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)

set(CMAKE_SKIP_BUILD_RPATH TRUE)   ###### <--------- This is absolutely necessary if you use linear algebra

find_package(openfpm 3.2.0 REQUIRED)
find_package(Threads)
find_package(MPI)

if (openfpm_FOUND)
	message("OpenFPM found")
	
	add_executable(sph_dlb  main.cpp)
	target_link_libraries(sph_dlb PUBLIC openfpm::binary_config)

	# or

	#target_include_directories(sph_dlb PUBLIC ${OPENFPM_INCLUDES})
	#target_compile_definitions(sph_dlb PUBLIC ${OPENFPM_DEFINITION})
	#target_link_libraries(sph_dlb PUBLIC ${OPENFPM_LIBS})
	#target_compile_options(sph_dlb PUBLIC ${OPENFPM_COMPILE_OPTIONS})
endif()

