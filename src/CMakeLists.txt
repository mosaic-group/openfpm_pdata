cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

########################### Executables


if(NOT CUDA_ON_BACKEND STREQUAL "None")
	set(CUDA_SOURCES 
	    Grid/tests/sgrid_dist_id_gpu_unit_tests.cu
	    Vector/cuda/vector_dist_gpu_MP_tests.cu
	    Vector/cuda/vector_dist_cuda_func_test.cu
	    Decomposition/cuda/decomposition_cuda_tests.cu
	    Vector/cuda/vector_dist_gpu_unit_tests.cu
	    Decomposition/cuda/Domain_icells_cart_unit_test.cu
		Amr/tests/amr_base_gpu_unit_tests.cu
		Grid/tests/grid_dist_id_unit_test.cu)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
        add_definitions("-DBOOST_MPL_CFG_HAS_TYPEOF")
endif()


if (CUDA_ON_BACKEND STREQUAL "OpenMP" OR CUDA_ON_BACKEND STREQUAL "SEQUENTIAL")
        set_source_files_properties(${CUDA_SOURCES} PROPERTIES COMPILE_FLAGS "-D__NVCC__ -DCUDART_VERSION=11000")
        set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE CXX)
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
                add_definitions("-x c++")
        endif()
endif()

if ( CUDA_ON_BACKEND STREQUAL "HIP" AND HIP_FOUND )

        list(APPEND HIP_HIPCC_FLAGS ${CMAKE_CXX_FLAGS_DEBUG})

        if (CMAKE_BUILD_TYPE STREQUAL "Debug")
                list(APPEND HIP_HIPCC_FLAGS -O0)
        endif()

	list(APPEND HIP_HIPCC_FLAGS -D__NVCC__ -D__HIP__  -DCUDART_VERSION=11000 -D__CUDACC__ -D__CUDACC_VER_MAJOR__=11 -D__CUDACC_VER_MINOR__=0 -D__CUDACC_VER_BUILD__=0 --std=c++14)
        set_source_files_properties(${CUDA_SOURCES} PROPERTIES LANGUAGE CXX)

	#        set(CMAKE_CXX_COMPILER ${HIP_HIPCC_EXECUTABLE})

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
							  )


    hip_add_library(ofpm_pdata STATIC lib/pdata.cpp)
    hip_add_library(ofpm_pdata_dl SHARED lib/pdata.cpp)

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
							  lib/pdata.cpp)

    add_library(ofpm_pdata STATIC lib/pdata.cpp)
    add_library(ofpm_pdata_dl SHARED lib/pdata.cpp)
	add_library(ofpm_pdata_python SHARED lib/pdata_python.cpp)

	set_property(TARGET pdata PROPERTY CUDA_ARCHITECTURES OFF)
    set_property(TARGET ofpm_pdata PROPERTY CUDA_ARCHITECTURES OFF)
    set_property(TARGET ofpm_pdata_dl PROPERTY CUDA_ARCHITECTURES OFF)
	set_property(TARGET ofpm_pdata_python PROPERTY CUDA_ARCHITECTURES OFF)

endif()

add_dependencies(pdata ofpmmemory)
add_dependencies(pdata vcluster)

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
#target_compile_options(pdata PRIVATE "-Wno-deprecated-declarations")
#target_compile_options(pdata PRIVATE "-Wno-undefined-var-template")
#target_compile_options(pdata PRIVATE "-Wno-macro-redefined")
#target_compile_options(pdata PRIVATE "-fsanitize=address")
#    target_link_options(pdata PRIVATE "-fsanitize=address")
    if (TEST_COVERAGE)
######################################## Add mem tests ############################

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
target_include_directories (pdata PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_include_directories (pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_devices/src/)
target_include_directories (pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_vcluster/src/)
target_include_directories (pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_data/src/)
target_include_directories (pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_io/src/)
target_include_directories (pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/config)
target_include_directories (pdata PUBLIC ${PETSC_INCLUDES})
target_include_directories (pdata PUBLIC ${LIBHILBERT_INCLUDE_DIRS})
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
if (OpenMP_FOUND)
	target_link_libraries(pdata OpenMP::OpenMP_CXX)
endif()

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

#add_definitions(-DMAKE_CELLLIST_DETERMINISTIC)

if (TEST_COVERAGE)
    target_link_libraries(pdata -lgcov --coverage)

    get_target_property(PDATA_SOURCES pdata SOURCES)
    get_target_property(MEM_SOURCES mem SOURCES)
    get_target_property(DATA_SOURCES mem_map SOURCES)
    get_target_property(VCLUSTER_SOURCES vcluster_test SOURCES)
    get_target_property(IO_SOURCES io SOURCES)
    get_target_property(NUMERIC_SOURCES numerics SOURCES)
    get_target_property(PDATA_INCLUDES pdata INCLUDE_DIRECTORIES)
    get_target_property(NUMERIC_INCLUDES numerics INCLUDE_DIRECTORIES)
    get_target_property(PDATA_LIBS pdata LINK_LIBRARIES)
    get_target_property(NUMERIC_LIBS numerics LINK_LIBRARIES)

    list(TRANSFORM MEM_SOURCES PREPEND "../openfpm_devices/src/")
    list(TRANSFORM DATA_SOURCES PREPEND "../openfpm_data/src/")
    list(TRANSFORM VCLUSTER_SOURCES PREPEND "../openfpm_vcluster/src/")
    list(TRANSFORM IO_SOURCES PREPEND "../openfpm_io/src/")
    list(TRANSFORM NUMERIC_SOURCES PREPEND "../openfpm_numerics/src/")

    set(PDATA_FULL_SOURCES ${PDATA_SOURCES} ${MEM_SOURCES} ${DATA_SOURCES} ${VCLUSTER_SOURCES} ${IO_SOURCES} ${NUMERIC_SOURCES})
    set(PDATA_FULL_SOURCES_CU ${PDATA_SOURCES} ${MEM_SOURCES} ${DATA_SOURCES} ${VCLUSTER_SOURCES} ${IO_SOURCES} ${NUMERIC_SOURCES})
    list(FILTER PDATA_FULL_SOURCES_CU INCLUDE REGEX ^.*\.cu)

    ### Test coverage are done in SEQUENTIAL
    if (NOT CUDA_ON_BACKEND STREQUAL "SEQUENTIAL")
	    message(FATAL_ERROR "Error test coverage must be done with CUDA_ON_BACKEND=SEQUENTIAL backend")
    endif()

    set_source_files_properties(${PDATA_FULL_SOURCES_CU} PROPERTIES COMPILE_FLAGS "-D__NVCC__ -DCUDART_VERSION=11000")
    set_source_files_properties(${PDATA_FULL_SOURCES_CU} PROPERTIES LANGUAGE CXX)
    add_definitions("-x c++")

    add_executable(pdata_full ${PDATA_FULL_SOURCES})
    
    target_compile_definitions(pdata_full PUBLIC NO_INIT_AND_MAIN)
    target_compile_options(pdata_full PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -fprofile-arcs -ftest-coverage>)
    target_link_options(pdata_full PUBLIC -fprofile-arcs -ftest-coverage)

    target_link_options(pdata_full PUBLIC -fprofile-arcs -ftest-coverage)

    if (ENABLE_ASAN)
        target_compile_options(pdata_full PUBLIC $<$<COMPILE_LANGUAGE:CUDA>: -Xcompiler "-fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g" >)
        target_compile_options(pdata_full PRIVATE $<$<COMPILE_LANGUAGE:CXX>: -fsanitize=address -fno-optimize-sibling-calls -fsanitize-address-use-after-scope -fno-omit-frame-pointer -g >)
    endif()

    target_link_libraries(pdata_full -lgcov --coverage)
    target_link_libraries(pdata_full ${PDATA_LIBS})
    target_link_libraries(pdata_full ${NUMERIC_LIBS})
    target_link_libraries(pdata_full Boost::thread)
    target_include_directories(pdata_full PUBLIC ${PDATA_INCLUDES})
    target_include_directories(pdata_full PUBLIC ${NUMERIC_INCLUDES})
endif()



target_include_directories (ofpm_pdata PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_include_directories (ofpm_pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (ofpm_pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_data/src/)
target_include_directories (ofpm_pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/config)
target_include_directories (ofpm_pdata PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_devices/src/)
target_include_directories (ofpm_pdata PUBLIC ${Boost_INCLUDE_DIRS})

target_include_directories (ofpm_pdata_dl PUBLIC ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_include_directories (ofpm_pdata_dl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (ofpm_pdata_dl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_data/src/)
target_include_directories (ofpm_pdata_dl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/config)
target_include_directories (ofpm_pdata_dl PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_devices/src/)
target_include_directories (ofpm_pdata_dl PUBLIC ${Boost_INCLUDE_DIRS})

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
install(TARGETS ofpm_pdata_dl EXPORT ofpm_pdata_dl_config  DESTINATION openfpm_pdata/lib COMPONENT OpenFPM)

########## Create openfpmConfig.cmake + openfpmConfigVersion.cmake

set(DOBC $ENV{DISABLE_OPENFPM_BINARY_CONFIG})

if(NOT DOBC EQUAL 1)

    add_library(binary_config INTERFACE)
    add_library(binary_config_ INTERFACE)
    add_library(binary_config_shared INTERFACE)
    add_library(binary_config_static INTERFACE)

    target_include_directories(
    binary_config_
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
    ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    ${PETSC_INCLUDES}
    ${HDF5_ROOT}/include
    ${LIBHILBERT_INCLUDE_DIRS}
    ${Vc_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
    )

    if (CUDA_ON_CPU)
        target_compile_definitions(binary_config_ INTERFACE CUDA_ON_CPU)
    endif()

    target_compile_options(binary_config_ INTERFACE $<$<COMPILE_LANGUAGE:CUDA>: ${WARNING_SUPPRESSION_AND_OPTION_NVCC} >)

    target_link_libraries(binary_config_ INTERFACE Boost::context)
    target_link_libraries(binary_config_ INTERFACE Boost::fiber)
    target_link_libraries(binary_config_ INTERFACE Boost::filesystem)
    target_link_libraries(binary_config_ INTERFACE Boost::iostreams)
    target_link_libraries(binary_config_ INTERFACE Boost::program_options)
    target_link_libraries(binary_config_ INTERFACE ${PARMETIS_LIBRARIES})
    target_link_libraries(binary_config_ INTERFACE -L${METIS_ROOT}/lib metis)
    target_link_libraries(binary_config_ INTERFACE ${HDF5_LIBRARIES})
    target_link_libraries(binary_config_ INTERFACE -L${LIBHILBERT_LIBRARY_DIRS} ${LIBHILBERT_LIBRARIES})
    target_link_libraries(binary_config_ INTERFACE ${PETSC_LIBRARIES})
    target_link_libraries(binary_config_ INTERFACE ${Vc_LIBRARIES})
    target_link_libraries(binary_config_ INTERFACE ${alpaka_LIBRARIES})
    target_link_libraries(binary_config_ INTERFACE ${MPI_C_LIBRARIES})

    # Not OK before CMake 3.13
    #target_link_libraries(binary_config INTERFACE $<INSTALL_PREFIX>/openfpm_vcluster/lib/$<TARGET_FILE_NAME:openfpm::vcluster> )
    #target_link_libraries(binary_config INTERFACE $<INSTALL_PREFIX>/openfpm_devices/lib/$<TARGET_FILE_NAME:openfpm::ofpmmemory> )
    target_link_libraries(binary_config_static INTERFACE ${CMAKE_INSTALL_PREFIX}/openfpm_vcluster/lib/libvcluster.a )
    target_link_libraries(binary_config_static INTERFACE ${CMAKE_INSTALL_PREFIX}/openfpm_devices/lib/libofpmmemory.a )
    target_link_libraries(binary_config_static INTERFACE ${CMAKE_INSTALL_PREFIX}/openfpm_pdata/lib/libofpm_pdata.a )

    target_link_libraries(binary_config_shared INTERFACE ${CMAKE_INSTALL_PREFIX}/openfpm_vcluster/lib/libvcluster_dl.so )
    target_link_libraries(binary_config_shared INTERFACE ${CMAKE_INSTALL_PREFIX}/openfpm_devices/lib/libofpmmemory_dl.so )
    target_link_libraries(binary_config_shared INTERFACE ${CMAKE_INSTALL_PREFIX}/openfpm_pdata/lib/libofpm_pdata_dl.so )

    target_link_libraries(binary_config_ INTERFACE ${CUDA_LIBRARIES} )

    # Does not work before Cmake 3.13
    #install(TARGETS binary_config vcluster ofpmmemory EXPORT openfpm_config  CONFIGURATIONS)
    install(TARGETS binary_config EXPORT openfpm_config  CONFIGURATIONS)
    install(TARGETS binary_config_ EXPORT openfpm_config  CONFIGURATIONS)
    install(TARGETS binary_config_static EXPORT openfpm_config  CONFIGURATIONS)
    install(TARGETS binary_config_shared EXPORT openfpm_config  CONFIGURATIONS)

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
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/openfpmConfig-Vars.cmake"
    DESTINATION
        cmake/
    )
    
endif()

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

