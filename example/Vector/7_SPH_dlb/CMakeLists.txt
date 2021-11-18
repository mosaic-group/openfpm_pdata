cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)

cmake_policy(SET CMP0074 NEW)

set(CMAKE_SKIP_BUILD_RPATH TRUE)   ###### <--------- This is absolutely necessary if you use linear algebra


find_package(Threads)
find_package(MPI)
find_package(Boost 1.75.0 COMPONENTS program_options iostreams filesystem fiber context)

# CMAKE_PREFIX_PATH can be used to control where to find OpenFPM, the other way is set the path with openfpm_ROOT and activate policy CMP0074
find_package(openfpm 3.2.0 REQUIRED)

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

