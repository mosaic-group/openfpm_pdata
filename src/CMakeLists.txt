cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

add_definitions(-DSCAN_WITH_CUB)

########################### Executables


if(CUDA_FOUND OR CUDA_ON_CPU OR HIP_FOUND)
	set(CUDA_SOURCES 
	    Grid/tests/sgrid_dist_id_gpu_unit_tests.cu
	    Vector/cuda/vector_dist_gpu_MP_tests.cu 
	    Vector/cuda/vector_dist_cuda_func_test.cu
	    Decomposition/cuda/decomposition_cuda_tests.cu
	    Vector/cuda/vector_dist_gpu_unit_tests.cu
	    Decomposition/cuda/Domain_icells_cart_unit_test.cu
	    Amr/tests/amr_base_gpu_unit_tests.cu)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        add_definitions("-DBOOST_MPL_CFG_HAS_TYPEOF")
endif()

if (CUDA_ON_CPU)
        add_definitions(-DCUDA_ON_CPU)
        set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE CXX)
	set_source_files_properties(${CUDA_SOURCES} PROPERTIES COMPILE_FLAGS "-D__NVCC__ -DCUDART_VERSION=11000")
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
                add_definitions("-x c++")
        endif()
endif()

if ( HIP_ENABLE AND HIP_FOUND )

        list(APPEND HIP_HIPCC_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})

        if (CMAKE_BUILD_TYPE STREQUAL "Debug")
                list(APPEND HIP_HIPCC_FLAGS -O0)
        endif()

	list(APPEND HIP_HIPCC_FLAGS -D__NVCC__ -D__HIP__  -DCUDART_VERSION=11000 -D__CUDACC__ -D__CUDACC_VER_MAJOR__=11 -D__CUDACC_VER_MINOR__=0 -D__CUDACC_VER_BUILD__=0)
        set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE CXX)

        set(CMAKE_CXX_COMPILER ${HIP_HIPCC_EXECUTABLE})

        hip_add_executable(pdata ${CUDA_SOURCES} ${OPENFPM_INIT_FILE} main.cpp
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
							  lib/pdata.cpp 
							  test_multiple_o.cpp
							  )
							  

	hip_add_library(ofpm_pdata STATIC lib/pdata.cpp)

else()

	add_executable(pdata ${OPENFPM_INIT_FILE} ${CUDA_SOURCES} main.cpp
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
							  lib/pdata.cpp test_multiple_o.cpp)

	add_library(ofpm_pdata STATIC lib/pdata.cpp)

endif()

add_dependencies(pdata ofpmmemory)
add_dependencies(pdata vcluster)

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(pdata PRIVATE "-Wno-undefined-var-template")

        # Instead of enabling this option let's check why is redefined
        #target_compile_options(pdata PRIVATE "-Wno-macro-redefined")
endif()
if ( CMAKE_COMPILER_IS_GNUCC )
        target_compile_options(pdata PRIVATE "-Wno-deprecated-declarations")
        target_compile_options(pdata PRIVATE "-Wno-undefined-var-template")

        # Instead of enabling this option let's check why is redefined
        #target_compile_options(pdata PRIVATE "-Wno-macro-redefined")
endif()

#add_executable(isolation_pdata ${OPENFPM_INIT_FILE} isolation.cu
#														  lib/pdata.cpp  
#														  ../openfpm_devices/src/memory/HeapMemory.cpp 
#														  ../openfpm_devices/src/memory/CudaMemory.cu
#														  ../openfpm_devices/src/memory/PtrMemory.cpp 
#														  ../openfpm_vcluster/src/VCluster/VCluster.cpp 
#														  )

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	target_compile_options(pdata PRIVATE "-Wno-undefined-var-template")
	target_compile_options(pdata PRIVATE "-Wno-macro-redefined")
endif()
if ( CMAKE_COMPILER_IS_GNUCC )
    target_compile_options(pdata PRIVATE "-Wno-deprecated-declarations")
	target_compile_options(pdata PRIVATE "-Wno-undefined-var-template")
	target_compile_options(pdata PRIVATE "-Wno-macro-redefined")
    #target_compile_options(pdata PRIVATE "-fsanitize=address")
    #    target_link_options(pdata PRIVATE "-fsanitize=address")
    if (TEST_COVERAGE)
        target_compile_options(pdata PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -fprofile-arcs -ftest-coverage>)
    endif()
endif()

if (CMAKE_CUDA_COMPILER_ID STREQUAL "Clang")
        add_definitions(-D__STRICT_ANSI__)
endif()

if (ENABLE_ASAN)
    target_compile_options(pdata PUBLIC $<$<COMPILE_LANGUAGE:CUDA>: -Xcompiler "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g" >)
    target_compile_options(pdata PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g >)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer")
	add_definitions(-DENABLE_ASAN)
endif()


add_test(NAME pdata_3_proc COMMAND mpirun -np 3 ./pdata)
add_test(NAME pdata_4_proc COMMAND mpirun -np 4 ./pdata)

###########################

if (CUDA_FOUND)
	target_include_directories (pdata PUBLIC ${MPI_C_INCLUDE_DIRS})
        if (TEST_COVERAGE)
		target_compile_options(pdata PRIVATE $<$<COMPILE_LANGUAGE:CUDA>: -Xcompiler "-fprofile-arcs -ftest-coverage" >)
        endif()
endif()

if(TEST_PERFORMANCE)
        target_include_directories (pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_numerics/src/)
endif()
target_include_directories (pdata PUBLIC ${PARMETIS_ROOT}/include)
target_include_directories (pdata PUBLIC ${METIS_ROOT}/include)
target_include_directories (pdata PUBLIC ${HDF5_ROOT}/include)
target_include_directories (pdata PUBLIC ${CUDA_INCLUDE_DIRS})
target_include_directories (pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_devices/src/)
target_include_directories (pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_vcluster/src/)
target_include_directories (pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_data/src/)
target_include_directories (pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_io/src/)
target_include_directories (pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/config)
target_include_directories (pdata PUBLIC ${PETSC_INCLUDES})
target_include_directories (pdata PUBLIC ${LIBHILBERT_INCLUDE_DIRS})
if(EIGEN3_FOUND)
	target_include_directories (pdata PUBLIC ${EIGEN3_INCLUDE_DIR})
endif()

target_include_directories (pdata PUBLIC ${ALPAKA_ROOT}/include)
target_include_directories (pdata PUBLIC ${Vc_INCLUDE_DIR})
target_include_directories (pdata PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories (pdata PUBLIC ${MPI_C_INCLUDE_DIRS})

#target_include_directories (isolation_pdata PUBLIC ${PARMETIS_ROOT}/include)
#target_include_directories (isolation_pdata PUBLIC ${METIS_ROOT}/include)
#target_include_directories (isolation_pdata PUBLIC ${CUDA_INCLUDE_DIRS})
#target_include_directories (isolation_pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
#target_include_directories (isolation_pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_devices/src/)
#target_include_directories (isolation_pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_vcluster/src/)
#target_include_directories (isolation_pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_data/src/)
#target_include_directories (isolation_pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_io/src/)
#target_include_directories (isolation_pdata PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/config)
#target_include_directories (isolation_pdata PUBLIC ${PETSC_INCLUDES})
#target_include_directories (isolation_pdata PUBLIC ${HDF5_ROOT}/include)
#target_include_directories (isolation_pdata PUBLIC ${LIBHILBERT_INCLUDE_DIRS})
#target_include_directories (isolation_pdata PUBLIC ${Vc_INCLUDE_DIR})
#target_include_directories (isolation_pdata PUBLIC ${Boost_INCLUDE_DIRS})


target_link_libraries(pdata ${Boost_LIBRARIES})
target_link_libraries(pdata ${PARMETIS_LIBRARIES})
target_link_libraries(pdata -L${METIS_ROOT}/lib metis)
target_link_libraries(pdata ${HDF5_LIBRARIES})
target_link_libraries(pdata -L${LIBHILBERT_LIBRARY_DIRS} ${LIBHILBERT_LIBRARIES})
target_link_libraries(pdata ${PETSC_LIBRARIES})
target_link_libraries(pdata ${Vc_LIBRARIES})
target_link_libraries(pdata ${alpaka_LIBRARIES})
target_link_libraries(pdata ${MPI_C_LIBRARIES})
target_link_libraries(pdata ${MPI_CXX_LIBRARIES})
target_link_libraries(pdata vcluster)
target_link_libraries(pdata ofpmmemory)

#target_link_libraries(isolation_pdata ${Boost_LIBRARIES})
#target_link_libraries(isolation_pdata ${PARMETIS_LIBRARIES})
#target_link_libraries(isolation_pdata -L${METIS_ROOT}/lib metis)
#target_link_libraries(isolation_pdata ${HDF5_LIBRARIES})
#target_link_libraries(isolation_pdata -L${LIBHILBERT_LIBRARY_DIRS} ${LIBHILBERT_LIBRARIES})
#target_link_libraries(isolation_pdata ${PETSC_LIBRARIES})
#target_link_libraries(isolation_pdata ${Vc_LIBRARIES})

if (TEST_PERFORMANCE)
	target_link_libraries(pdata  ${Boost_FILESYSTEM_LIBRARY})
        target_link_libraries(pdata ${Boost_SYSTEM_LIBRARY})
endif()

add_definitions(-DSCAN_WITH_CUB)
#add_definitions(-DMAKE_CELLLIST_DETERMINISTIC)

if (TEST_COVERAGE)
    target_link_libraries(pdata -lgcov --coverage)
endif()

# Debug!
# Hack found at https://github.com/LLNL/scr/issues/130#issuecomment-402815952
IF(MPI_CXX_FOUND)
	INCLUDE_DIRECTORIES(${MPI_CXX_INCLUDE_PATH})
#	LIST(APPEND SCR_EXTERNAL_LIBS ${MPI_CXX_LIBRARIES})
	target_link_libraries(pdata ${MPI_CXX_LIBRARIES})
ENDIF(MPI_CXX_FOUND)

target_include_directories (ofpm_pdata PUBLIC ${CUDA_INCLUDE_DIRS})
target_include_directories (ofpm_pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (ofpm_pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_data/src/)
target_include_directories (ofpm_pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/config)
target_include_directories (ofpm_pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_devices/src/)
target_include_directories (ofpm_pdata PUBLIC ${Boost_INCLUDE_DIRS})

target_compile_definitions(pdata PRIVATE ${MPI_VENDOR})

if(PETSC_FOUND)
	target_link_libraries(pdata ${PETSC_LIBRARIES})
endif()



# Request that particles be built with -std=c++11
# As this is a public compile feature anything that links to particles
# will also build with -std=c++11
target_compile_features(pdata PUBLIC cxx_std_14)
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
	      DESTINATION openfpm_pdata/include/Decomposition/ 
	      COMPONENT OpenFPM)

install(FILES Decomposition/Distribution/metis_util.hpp 
	      Decomposition/Distribution/SpaceDistribution.hpp 
	      Decomposition/Distribution/parmetis_dist_util.hpp  
	      Decomposition/Distribution/parmetis_util.hpp 
	      Decomposition/Distribution/MetisDistribution.hpp 
	      Decomposition/Distribution/ParMetisDistribution.hpp 
	      Decomposition/Distribution/DistParMetisDistribution.hpp
	      Decomposition/Distribution/BoxDistribution.hpp  
	      DESTINATION openfpm_pdata/include/Decomposition/Distribution 
	      COMPONENT OpenFPM)

install(FILES Decomposition/cuda/ie_ghost_gpu.cuh
	      Decomposition/cuda/CartDecomposition_gpu.cuh
	      DESTINATION openfpm_pdata/include/Decomposition/cuda 
	      COMPONENT OpenFPM)


install(FILES Grid/grid_dist_id.hpp 
	      Grid/grid_dist_id_comm.hpp
	      Grid/grid_dist_util.hpp  
	      Grid/grid_dist_key.hpp 
	      Grid/staggered_dist_grid.hpp 
	      Grid/staggered_dist_grid_util.hpp 
	      Grid/staggered_dist_grid_copy.hpp
	      DESTINATION openfpm_pdata/include/Grid/ 
	      COMPONENT OpenFPM)

install(FILES Grid/cuda/grid_dist_id_kernels.cuh
	      Grid/cuda/grid_dist_id_iterator_gpu.cuh
	DESTINATION openfpm_pdata/include/Grid/cuda/ 
	COMPONENT OpenFPM)

install(FILES Amr/grid_dist_amr_key_iterator.hpp 
	      Amr/grid_dist_amr_key.hpp
	      Amr/grid_dist_amr.hpp
	      DESTINATION openfpm_pdata/include/Amr/ 
	      COMPONENT OpenFPM)

install(FILES Grid/Iterators/grid_dist_id_iterator_util.hpp
              Grid/Iterators/grid_dist_id_iterator_dec.hpp
              Grid/Iterators/grid_dist_id_iterator_dec_skin.hpp
              Grid/Iterators/grid_dist_id_iterator_sub.hpp
	      Grid/Iterators/grid_dist_id_iterator.hpp
	      DESTINATION openfpm_pdata/include/Grid/Iterators 
	      COMPONENT OpenFPM)


install(FILES Vector/se_class3_vector.hpp  
	      Vector/vector_dist_multiphase_functions.hpp 
	      Vector/vector_dist_comm.hpp Vector/vector_dist.hpp 
	      Vector/vector_dist_ofb.hpp 
	      Vector/vector_dist_key.hpp
	      Vector/vector_dist_kernel.hpp
	      Vector/vector_dist_subset.hpp
	      DESTINATION openfpm_pdata/include/Vector
	      COMPONENT OpenFPM)

install(FILES util/common_pdata.hpp
	      DESTINATION openfpm_pdata/include/util
	      COMPONENT OpenFPM)

install(FILES Vector/Iterators/vector_dist_iterator.hpp
	      DESTINATION openfpm_pdata/include/Vector/Iterators/ 
	      COMPONENT OpenFPM)

install(FILES Vector/util/vector_dist_funcs.hpp
	      DESTINATION openfpm_pdata/include/Vector/util 
	      COMPONENT OpenFPM)

install(FILES Vector/cuda/vector_dist_comm_util_funcs.cuh
	      Vector/cuda/vector_dist_cuda_funcs.cuh
	      Vector/cuda/vector_dist_operators_list_ker.hpp
	DESTINATION openfpm_pdata/include/Vector/cuda 
	COMPONENT OpenFPM)

install(FILES Graph/ids.hpp Graph/dist_map_graph.hpp 
	      Graph/DistGraphFactory.hpp
              DESTINATION openfpm_pdata/include/Graph 
	      COMPONENT OpenFPM)

install(FILES example.mk
	      SubdomainGraphNodes.hpp
              DESTINATION openfpm_pdata/include/ )

install(FILES DLB/DLB.hpp DLB/LB_Model.hpp
	DESTINATION openfpm_pdata/include/DLB 
	COMPONENT OpenFPM)

install(FILES config/config.h
        DESTINATION openfpm_pdata/include/config 
	COMPONENT OpenFPM)

install(FILES lib/pdata.hpp
        DESTINATION openfpm_pdata/include/lib 
	COMPONENT OpenFPM)

install(FILES Debug/debug.hpp
	DESTINATION openfpm_pdata/include/Debug 
	COMPONENT OpenFPM)

install(TARGETS ofpm_pdata EXPORT ofpm_pdata_config  DESTINATION openfpm_pdata/lib COMPONENT OpenFPM)

########## Create openfpmConfig.cmake + openfpmConfigVersion.cmake

add_library(binary_config INTERFACE)

target_include_directories(
  binary_config
  INTERFACE 
  ${CMAKE_INSTALL_PREFIX}/openfpm_pdata/include
  ${CMAKE_INSTALL_PREFIX}/openfpm_data/include
  ${CMAKE_INSTALL_PREFIX}/openfpm_pdata/include/config
  ${CMAKE_INSTALL_PREFIX}/openfpm_io/include
  ${CMAKE_INSTALL_PREFIX}/openfpm_vcluster/include
  ${CMAKE_INSTALL_PREFIX}/openfpm_devices/include
  ${CMAKE_INSTALL_PREFIX}/openfpm_numerics/include
  ${PARMETIS_ROOT}/include
  ${METIS_ROOT}/include
  ${CUDA_INCLUDE_DIRS}
  ${PETSC_INCLUDES}
  ${HDF5_ROOT}/include
  ${LIBHILBERT_INCLUDE_DIRS}
  ${Vc_INCLUDE_DIR}
  ${Boost_INCLUDE_DIRS}
  )

if (CUDA_ON_CPU)
	target_compile_definitions(binary_config INTERFACE CUDA_ON_CPU)
endif()

target_compile_options(binary_config INTERFACE $<$<COMPILE_LANGUAGE:CUDA>: ${WARNING_SUPPRESSION_AND_OPTION_NVCC} >)

target_link_libraries(binary_config INTERFACE ${Boost_LIBRARIES})
target_link_libraries(binary_config INTERFACE ${PARMETIS_LIBRARIES})
target_link_libraries(binary_config INTERFACE -L${METIS_ROOT}/lib metis)
target_link_libraries(binary_config INTERFACE ${HDF5_LIBRARIES})
target_link_libraries(binary_config INTERFACE -L${LIBHILBERT_LIBRARY_DIRS} ${LIBHILBERT_LIBRARIES})
target_link_libraries(binary_config INTERFACE ${PETSC_LIBRARIES})
target_link_libraries(binary_config INTERFACE ${Vc_LIBRARIES})
target_link_libraries(binary_config INTERFACE ${alpaka_LIBRARIES})
target_link_libraries(binary_config INTERFACE ${MPI_C_LIBRARIES})

# Not OK before CMake 3.13
#target_link_libraries(binary_config INTERFACE $<INSTALL_PREFIX>/openfpm_vcluster/lib/$<TARGET_FILE_NAME:openfpm::vcluster> )
#target_link_libraries(binary_config INTERFACE $<INSTALL_PREFIX>/openfpm_devices/lib/$<TARGET_FILE_NAME:openfpm::ofpmmemory> )
target_link_libraries(binary_config INTERFACE ${CMAKE_INSTALL_PREFIX}/openfpm_vcluster/lib/libvcluster.a )
target_link_libraries(binary_config INTERFACE ${CMAKE_INSTALL_PREFIX}/openfpm_devices/lib/libofpmmemory.a )
target_link_libraries(binary_config INTERFACE ${CMAKE_INSTALL_PREFIX}/openfpm_pdata/lib/libofpm_pdata.a )
target_link_libraries(binary_config INTERFACE ${CUDA_LIBRARIES} )

# Does not work before Cmake 3.13
#install(TARGETS binary_config vcluster ofpmmemory EXPORT openfpm_config  CONFIGURATIONS)
install(TARGETS binary_config EXPORT openfpm_config  CONFIGURATIONS)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/openfpm_cmake/openfpmConfigVersion.cmake"
  VERSION ${openfpm_VERSION}
  COMPATIBILITY AnyNewerVersion
)


install(EXPORT openfpm_config
        DESTINATION cmake/
        NAMESPACE openfpm::
        FILE openfpmConfig.cmake)

install(
  FILES
    "${CMAKE_CURRENT_BINARY_DIR}/openfpm_cmake/openfpmConfigVersion.cmake"
  DESTINATION
    cmake/
)

install(
  FILES
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/openfpmConfigVars-configure.cmake"
  DESTINATION
    cmake/
)

#####################################################################

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

