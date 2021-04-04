cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(openfpm_pdata LANGUAGES C CXX)

find_package(openfpm 3.2.0 REQUIRED)

enable_language(CUDA)

if (openfpm_FOUND)
	message("OpenFPM found")
	
	add_executable(sph_dlb  main.cu)
	target_include_directories(sph_dlb PUBLIC ${OPENFPM_INCLUDES})
	target_compile_definitions(sph_dlb PUBLIC ${OPENFPM_DEFINITION})
	target_link_libraries(sph_dlb PUBLIC ${OPENFPM_LIBS})
	target_compile_options(sph_dlb PUBLIC ${OPENFPM_COMPILE_OPTIONS})
endif()

