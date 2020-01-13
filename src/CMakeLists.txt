cmake_minimum_required(VERSION 3.8 FATAL_ERROR)


########################### Executables


if ( HIP_ENABLE OR CUDA_FOUND)

	set(CUDA_SOURCES 
		Grid/tests/sgrid_dist_id_gpu_unit_tests.cu
		Vector/cuda/vector_dist_gpu_MP_tests.cu 
		Vector/cuda/vector_dist_cuda_func_test.cu
		Decomposition/cuda/decomposition_cuda_tests.cu
		Vector/cuda/vector_dist_gpu_unit_tests.cu
		../openfpm_devices/src/memory/CudaMemory.cu
		Decomposition/cuda/Domain_icells_cart_unit_test.cu
		Amr/tests/amr_base_gpu_unit_tests.cu)

else()
	set(CUDA_SOURCES)
endif()

set(CPP_SOURCES main.cpp
                                                                                                                  Amr/grid_dist_amr_unit_tests.cpp
                                                                                                                  Amr/tests/amr_base_unit_tests.cpp
                                                                                                                  Debug/debug_test.cpp
                                                                                                                  Grid/tests/grid_dist_id_HDF5_chckpnt_restart_test.cpp
                                                                                                                  Grid/tests/grid_dist_id_unit_test.cpp
                                                                                                                  Grid/tests/sgrid_dist_id_unit_tests.cpp
                                                                                                                  Grid/tests/grid_dist_id_dlb_unit_test.cpp
                                                                                                                  Grid/tests/staggered_grid_dist_unit_test.cpp
                                                                                                                  Vector/tests/vector_dist_cell_list_tests.cpp
                                                                                                                  Vector/tests/vector_dist_complex_prp_unit_test.cpp
                                                                                                                  Vector/tests/vector_dist_HDF5_chckpnt_restart_test.cpp
                                                                                                                  Vector/tests/vector_dist_MP_unit_tests.cpp
                                                                                                                  Vector/tests/vector_dist_NN_tests.cpp
                                                                                                                  Vector/tests/vector_dist_unit_test.cpp
                                                                                                                  pdata_performance.cpp
                                                                                                                  Decomposition/tests/CartDecomposition_unit_test.cpp
                                                                                                                  Decomposition/tests/shift_vect_converter_tests.cpp
                                                                                                                  Vector/performance/vector_dist_performance_util.cpp
                                                                                                                  lib/pdata.cpp test_multiple_o.cpp
                                                                                                                  ../openfpm_devices/src/memory/HeapMemory.cpp
                                                                                                                  ../openfpm_devices/src/memory/PtrMemory.cpp
                                                                                                                  ../openfpm_vcluster/src/VCluster/VCluster.cpp
                                                                                                                  ../openfpm_devices/src/Memleak_check.cpp)

if ( HIP_ENABLE AND HIP_FOUND )

	set(CMAKE_CXX_COMPILER ${HIP_HIPCC_EXECUTABLE})


	hip_add_executable(pdata ${OPENFPM_INIT_FILE} ${CUDA_SOURCES} ${CPP_SOURCES})
	hip_add_library(ofpm_pdata STATIC lib/pdata.cpp)

	target_include_directories (pdata PUBLIC ${HIP_ROOT_DIR}/include)
	#target_include_directories (pdata PUBLIC ${MPI_C_INCLUDE_DIRS}) # Well ........ it does not work, work around line below
	set_property(TARGET pdata PROPERTY COMPILE_FLAGS "-I${MPI_C_INCLUDE_DIRS}")
else()

	add_executable(pdata ${OPENFPM_INIT_FILE} ${CUDA_SOURCES} ${CPP_SOURCES})
	add_library(ofpm_pdata STATIC lib/pdata.cpp)

endif()

if ( CMAKE_COMPILER_IS_GNUCC )
    target_compile_options(pdata PRIVATE "-Wno-deprecated-declarations")
    if (TEST_COVERAGE)
        target_compile_options(pdata PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -fprofile-arcs -ftest-coverage>)
    endif()
endif()



add_test(NAME pdata_3_proc COMMAND mpirun -np 3 ./pdata)
add_test(NAME pdata_4_proc COMMAND mpirun -np 4 ./pdata)

###########################

if (CUDA_FOUND)
	target_compile_options(pdata PUBLIC $<$<COMPILE_LANGUAGE:CUDA>: ${WARNING_SUPPRESSION_AND_OPTION_NVCC} >)
	target_include_directories (pdata PUBLIC ${MPI_C_INCLUDE_DIRS})
        if (TEST_COVERAGE)
		target_compile_options(pdata PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -Xcompiler "-fprofile-arcs -ftest-coverage">)
        endif()
endif()


if(TEST_PERFORMANCE)
        target_include_directories (pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_numerics/src/)
endif()
target_include_directories (pdata PUBLIC ${PARMETIS_ROOT}/include)
target_include_directories (pdata PUBLIC ${METIS_ROOT}/include)
target_include_directories (pdata PUBLIC ${CUDA_INCLUDE_DIRS})
target_include_directories (pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_devices/src/)
target_include_directories (pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_vcluster/src/)
target_include_directories (pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_data/src/)
target_include_directories (pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_io/src/)
target_include_directories (pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/config)
target_include_directories (pdata PUBLIC ${PETSC_INCLUDES})
target_include_directories (pdata PUBLIC ${HDF5_ROOT}/include)
target_include_directories (pdata PUBLIC ${LIBHILBERT_INCLUDE_DIRS})
target_include_directories (pdata PUBLIC ${Boost_INCLUDE_DIRS})

if (HIP_ENABLE)
        SET(CMAKE_EXE_LINKER_FLAGS  "--amdgpu-target=${AMD_ARCH_COMPILE}")
endif()

target_link_libraries(pdata ${Boost_LIBRARIES})
target_link_libraries(pdata ${PARMETIS_LIBRARIES})
target_link_libraries(pdata -L${METIS_ROOT}/lib metis)
target_link_libraries(pdata -L${HDF5_ROOT}/lib hdf5 hdf5_hl)
target_link_libraries(pdata -L${LIBHILBERT_LIBRARY_DIRS} ${LIBHILBERT_LIBRARIES})
target_link_libraries(pdata ${PETSC_LIBRARIES})

if (TEST_PERFORMANCE)
	target_link_libraries(pdata  ${Boost_FILESYSTEM_LIBRARY})
        target_link_libraries(pdata ${Boost_SYSTEM_LIBRARY})
endif()

if (TEST_COVERAGE)
    target_link_libraries(pdata -lgcov --coverage)
endif()



target_include_directories (ofpm_pdata PUBLIC ${CUDA_INCLUDE_DIRS})
target_include_directories (ofpm_pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (ofpm_pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_data/src/)
target_include_directories (ofpm_pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/config)
target_include_directories (ofpm_pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_devices/src/)
target_include_directories (ofpm_pdata PUBLIC ${Boost_INCLUDE_DIRS})

target_compile_definitions(pdata PRIVATE ${MPI_VENDOR})

if(PETSC_FOUND)
	target_link_libraries(pdata ${PETSC_LIBRARIES})
endif()


# Request that particles be built with -std=c++11
# As this is a public compile feature anything that links to particles
# will also build with -std=c++11
target_compile_features(pdata PUBLIC cxx_std_11)
target_link_libraries(pdata ${MPI_C_LIBRARIES})
target_link_libraries(pdata m)
target_link_libraries(pdata c)
if (NOT APPLE)
    target_link_libraries(pdata rt)
endif ()

install(FILES Decomposition/CartDecomposition.hpp
       	      Decomposition/Domain_icells_cart.hpp	
	      Decomposition/shift_vect_converter.hpp 
	      Decomposition/CartDecomposition_ext.hpp  
	      Decomposition/common.hpp 
	      Decomposition/Decomposition.hpp  
	      Decomposition/ie_ghost.hpp
              Decomposition/Domain_NN_calculator_cart.hpp 
	      Decomposition/nn_processor.hpp Decomposition/ie_loc_ghost.hpp 
	      Decomposition/ORB.hpp
	      Decomposition/dec_optimizer.hpp
	      DESTINATION openfpm_pdata/include/Decomposition/ )

install(FILES Decomposition/Distribution/metis_util.hpp 
	      Decomposition/Distribution/SpaceDistribution.hpp 
	      Decomposition/Distribution/parmetis_dist_util.hpp  
	      Decomposition/Distribution/parmetis_util.hpp 
	      Decomposition/Distribution/MetisDistribution.hpp 
	      Decomposition/Distribution/ParMetisDistribution.hpp 
	      Decomposition/Distribution/DistParMetisDistribution.hpp  
	      DESTINATION openfpm_pdata/include/Decomposition/Distribution )

install(FILES Decomposition/cuda/ie_ghost_gpu.cuh
	      Decomposition/cuda/CartDecomposition_gpu.cuh
	      DESTINATION openfpm_pdata/include/Decomposition/cuda )

install(FILES Grid/grid_dist_id.hpp 
	      Grid/grid_dist_id_comm.hpp
	      Grid/grid_dist_util.hpp  
	      Grid/grid_dist_key.hpp 
	      Grid/staggered_dist_grid.hpp 
	      Grid/staggered_dist_grid_util.hpp 
	      Grid/staggered_dist_grid_copy.hpp
	      DESTINATION openfpm_pdata/include/Grid/ )

install(FILES Grid/Iterators/grid_dist_id_iterator_util.hpp
              Grid/Iterators/grid_dist_id_iterator_dec.hpp
              Grid/Iterators/grid_dist_id_iterator_dec_skin.hpp
              Grid/Iterators/grid_dist_id_iterator_sub.hpp
	      Grid/Iterators/grid_dist_id_iterator.hpp
	      DESTINATION openfpm_pdata/include/Grid/Iterators )


install(FILES Vector/se_class3_vector.hpp  
	      Vector/vector_dist_multiphase_functions.hpp 
	      Vector/vector_dist_comm.hpp Vector/vector_dist.hpp 
	      Vector/vector_dist_ofb.hpp 
	      Vector/vector_dist_key.hpp
	      Vector/vector_dist_kernel.hpp
	      DESTINATION openfpm_pdata/include/Vector )

install(FILES Vector/Iterators/vector_dist_iterator.hpp
	      DESTINATION openfpm_pdata/include/Vector/Iterators/ )

install(FILES Vector/util/vector_dist_funcs.hpp
	      DESTINATION openfpm_pdata/include/Vector/util )

install(FILES Vector/cuda/vector_dist_comm_util_funcs.cuh
	      Vector/cuda/vector_dist_cuda_funcs.cuh
	      Vector/cuda/vector_dist_operators_list_ker.hpp
	DESTINATION openfpm_pdata/include/Vector/cuda )

install(FILES Graph/ids.hpp Graph/dist_map_graph.hpp 
	      Graph/DistGraphFactory.hpp
              DESTINATION openfpm_pdata/include/Graph )

install(FILES example.mk
	      SubdomainGraphNodes.hpp
              DESTINATION openfpm_pdata/include/ )

install(FILES DLB/DLB.hpp DLB/LB_Model.hpp
	DESTINATION openfpm_pdata/include/DLB )

install(FILES config/config.h
        DESTINATION openfpm_pdata/include/config )

install(FILES lib/pdata.hpp
        DESTINATION openfpm_pdata/include/lib )

install(FILES Debug/debug.hpp
	DESTINATION openfpm_pdata/include/Debug )

install(TARGETS ofpm_pdata DESTINATION openfpm_pdata/lib)

#if(BUILD_TESTING)

#  add_executable(particle_test test.cu)

#  set_target_properties(particle_test PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#  target_link_libraries(particle_test PRIVATE particles)

#  add_test(NAME particles_10k COMMAND particle_test 10000 )
#  add_test(NAME particles_256k COMMAND particle_test 256000 )

#  if(APPLE)
    # We need to add the default path to the driver (libcuda.dylib) as an rpath,
    # so that the static cuda runtime can find it at runtime.
    #    set_property(TARGET particle_test PROPERTY BUILD_RPATH ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
    #  endif()
    #endif()

