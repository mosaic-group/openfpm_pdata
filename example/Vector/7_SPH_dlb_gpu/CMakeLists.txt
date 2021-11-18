cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)

cmake_policy(SET CMP0074 NEW)

set(CMAKE_SKIP_BUILD_RPATH TRUE)   ###### <--------- This is absolutely necessary if you use linear algebra

find_package(Threads)
find_package(MPI)
find_package(Boost 1.75.0 COMPONENTS program_options iostreams filesystem fiber context)
# CMAKE_PREFIX_PATH can be used to control where to find OpenFPM, the other way is set the path with openfpm_ROOT and activate policy CMP0074
find_package(openfpm 3.2.0 REQUIRED)
enable_language(CUDA)

if (openfpm_FOUND)
	
	# if we are not using CUDA as backend, but SEQUENTIAL or OpenMP you habe to uncomment this lines
#	set_source_files_properties(main.cu PROPERTIES COMPILE_FLAGS "-D__NVCC__ -DCUDART_VERSION=11000")  # <---- ENABLE IF YOU USE CUDA ON BACKEND SEQUENTIAL or OpenMP
#	set_source_files_properties(main.cu PROPERTIES LANGUAGE CXX)                                       # <---- ENABLE IF YOU USE CUDA ON BACKEND SEQUENTIAL or OpenMP
#	add_definitions("-x c++")                                                                          # <---- ENABLE IF YOU USE CUDA ON BACKEND SEQUENTIAL or OpenMP

	add_executable(sph_dlb  main.cu)
	target_link_libraries(sph_dlb PUBLIC openfpm::binary_config)
	set_property(TARGET sph_dlb PROPERTY CUDA_ARCHITECTURES OFF)	

	#or

	#target_include_directories(sph_dlb PUBLIC ${OPENFPM_INCLUDES})
	#target_compile_definitions(sph_dlb PUBLIC ${OPENFPM_DEFINITION})
	#target_link_libraries(sph_dlb PUBLIC ${OPENFPM_LIBS})
	#target_compile_options(sph_dlb PUBLIC ${OPENFPM_COMPILE_OPTIONS})
endif()

