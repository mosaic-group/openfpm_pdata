cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

########################### Executables

if(CUDA_FOUND)
	set(CUDA_SOURCES ../openfpm_devices/src/memory/CudaMemory.cu)
endif()

if ( HIP_ENABLE AND HIP_FOUND )

        set(CMAKE_CXX_COMPILER ${HIP_HIPCC_EXECUTABLE})

        hip_add_executable(cart_dec CartDecomposition_gen_vtk.cpp ../src/lib/pdata.cpp ${CUDA_SOURCES} ../openfpm_devices/src/memory/HeapMemory.cpp  ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp)

        hip_add_executable(metis_dec Metis_gen_vtk.cpp ../src/lib/pdata.cpp ${CUDA_SOURCES}  ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp)

        hip_add_executable(dom_box  domain_gen_vtk.cpp ../src/lib/pdata.cpp ${CUDA_SOURCES}  ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp)

        hip_add_executable(vector_dist vector.cpp ${CUDA_SOURCES}  ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_vcluster/src/VCluster/VCluster.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_devices/src/Memleak_check.cpp)

        #target_include_directories (pdata PUBLIC ${MPI_C_INCLUDE_DIRS}) # Well ........ it does not work, work around line below
        set_property(TARGET cart_dec PROPERTY COMPILE_FLAGS "-I${MPI_C_INCLUDE_DIRS}")
	set_property(TARGET metis_dec PROPERTY COMPILE_FLAGS "-I${MPI_C_INCLUDE_DIRS}")
	set_property(TARGET dom_box PROPERTY COMPILE_FLAGS "-I${MPI_C_INCLUDE_DIRS}")
	set_property(TARGET vector_dist PROPERTY COMPILE_FLAGS "-I${MPI_C_INCLUDE_DIRS}")

else()

	add_executable(cart_dec CartDecomposition_gen_vtk.cpp ../src/lib/pdata.cpp ${CUDA_SOURCES} ../openfpm_devices/src/memory/HeapMemory.cpp  ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp)

	add_executable(metis_dec Metis_gen_vtk.cpp ../src/lib/pdata.cpp ${CUDA_SOURCES}  ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp)

	add_executable(dom_box  domain_gen_vtk.cpp ../src/lib/pdata.cpp ${CUDA_SOURCES}  ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp)

	add_executable(vector_dist vector.cpp ${CUDA_SOURCES}  ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_vcluster/src/VCluster/VCluster.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_devices/src/Memleak_check.cpp)

endif()

if (HIP_FOUND)
        SET(CMAKE_EXE_LINKER_FLAGS  "--amdgpu-target=${AMD_ARCH_COMPILE}")
endif()

###########################


include_directories (${CUDA_INCLUDE_DIRS})
include_directories (${CMAKE_CURRENT_SOURCE_DIR})
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/../src/)
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_devices/src/)
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_vcluster/src/)
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_data/src/)
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/../openfpm_io/src/)
include_directories (${CMAKE_CURRENT_SOURCE_DIR}/../src/config)
include_directories (${HIP_ROOT_DIR}/include)
if (PETSC_FOUND)
	include_directories (${PETSC_INCLUDES})
endif()
include_directories (${PARMETIS_ROOT}/include)
include_directories (${METIS_ROOT}/include)
include_directories (${Boost_INCLUDE_DIRS})

#include_directories (${HDF5_INCLUDE_DIRS})
include_directories(${HDF5_ROOT}/include)
include_directories (${LIBHILBERT_INCLUDE_DIRS})

target_link_libraries(cart_dec ${Boost_LIBRARIES})
target_link_libraries(cart_dec -L${PARMETIS_ROOT}/lib parmetis)
target_link_libraries(cart_dec -L${METIS_ROOT}/lib metis)
target_link_libraries(cart_dec -L${HDF5_ROOT}/lib hdf5 hdf5_hl)
#target_link_libraries(cart_dec ${HDF5_LIBRARIES})
target_link_libraries(cart_dec -L${LIBHILBERT_LIBRARY_DIRS} ${LIBHILBERT_LIBRARIES})
target_link_libraries(cart_dec ${MPI_C_LIBRARIES})

target_link_libraries(metis_dec ${Boost_LIBRARIES})
target_link_libraries(metis_dec -L${PARMETIS_ROOT}/lib parmetis)
target_link_libraries(metis_dec -L${METIS_ROOT}/lib metis)
target_link_libraries(metis_dec -L${HDF5_ROOT}/lib hdf5 hdf5_hl)
#target_link_libraries(metis_dec ${HDF5_LIBRARIES})
target_link_libraries(metis_dec -L${LIBHILBERT_LIBRARY_DIRS} ${LIBHILBERT_LIBRARIES})
target_link_libraries(metis_dec ${MPI_C_LIBRARIES})

target_link_libraries(dom_box ${Boost_LIBRARIES})
target_link_libraries(dom_box -L${PARMETIS_ROOT}/lib parmetis)
target_link_libraries(dom_box -L${METIS_ROOT}/lib metis)
target_link_libraries(dom_box -L${HDF5_ROOT}/lib hdf5 hdf5_hl)
#target_link_libraries(dom_box ${HDF5_LIBRARIES})
target_link_libraries(dom_box -L${LIBHILBERT_LIBRARY_DIRS} ${LIBHILBERT_LIBRARIES})
target_link_libraries(dom_box ${MPI_C_LIBRARIES})

target_link_libraries(vector_dist ${Boost_LIBRARIES})
target_link_libraries(vector_dist -L${PARMETIS_ROOT}/lib parmetis)
target_link_libraries(vector_dist -L${METIS_ROOT}/lib metis)
target_link_libraries(vector_dist -L${HDF5_ROOT}/lib hdf5 hdf5_hl)
#target_link_libraries(vector_dist ${HDF5_LIBRARIES})
target_link_libraries(vector_dist -L${LIBHILBERT_LIBRARY_DIRS} ${LIBHILBERT_LIBRARIES})
target_link_libraries(vector_dist ${MPI_C_LIBRARIES})

if(PETSC_FOUND)
	target_link_libraries(cart_dec  ${PETSC_LIBRARIES})
	target_link_libraries(metis_dec ${PETSC_LIBRARIES})
	target_link_libraries(dom_box ${PETSC_LIBRARIES})
	target_link_libraries(vector_dist ${PETSC_LIBRARIES})
endif()


# Request that particles be built with -std=c++11
# As this is a public compile feature anything that links to particles
# will also build with -std=c++11
target_compile_features(cart_dec PUBLIC cxx_std_11)
target_compile_features(metis_dec PUBLIC cxx_std_11)
target_compile_features(dom_box PUBLIC cxx_std_11)
target_compile_features(vector_dist PUBLIC cxx_std_11)

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

