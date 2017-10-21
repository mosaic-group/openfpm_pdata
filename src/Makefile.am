LINKLIBS = $(HDF5_LDFLAGS)  $(HDF5_LIBS) $(OPENMP_LDFLAGS) $(LIBHILBERT_LIB)  $(METIS_LIB) $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_IOSTREAMS_LIB) $(CUDA_LIBS) $(PETSC_LIB) $(PARMETIS_LIB) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) $(BOOST_CHRONO_LIB) $(BOOST_TIMER_LIB) $(BOOST_SYSTEM_LIB) $(LIBIFCORE)

noinst_PROGRAMS = pdata
pdata_SOURCES = main.cpp pdata_performance.cpp Grid/tests/grid_dist_id_unit_test.cpp Amr/grid_dist_amr_unit_tests.cpp lib/pdata.cpp  Amr/tests/amr_base_unit_tests.cpp  ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp
pdata_CXXFLAGS = $(HDF5_CPPFLAGS) $(OPENMP_CFLAGS) $(AM_CXXFLAGS) $(LIBHILBERT_INCLUDE) $(PETSC_INCLUDE) $(CUDA_CFLAGS) $(INCLUDES_PATH) $(PARMETIS_INCLUDE) $(METIS_INCLUDE) $(BOOST_CPPFLAGS) $(H5PART_INCLUDE) -DPARALLEL_IO  -Wno-unused-local-typedefs
pdata_CFLAGS = $(CUDA_CFLAGS)
pdata_LDADD = $(LINKLIBS) -lparmetis -lmetis
nobase_include_HEADERS = Decomposition/CartDecomposition.hpp Decomposition/CartDecomposition_ext.hpp  Decomposition/common.hpp Decomposition/Decomposition.hpp  Decomposition/ie_ghost.hpp \
         Decomposition/Domain_NN_calculator_cart.hpp Decomposition/nn_processor.hpp Decomposition/ie_loc_ghost.hpp Decomposition/ORB.hpp \
         Graph/CartesianGraphFactory.hpp \
         Grid/grid_dist_id.hpp Grid/grid_dist_id_comm.hpp Grid/Iterators/grid_dist_id_iterator_util.hpp Grid/Iterators/grid_dist_id_iterator_dec.hpp Grid/Iterators/grid_dist_id_iterator_dec_skin.hpp Grid/grid_dist_util.hpp  Grid/Iterators/grid_dist_id_iterator_sub.hpp Grid/Iterators/grid_dist_id_iterator.hpp Grid/grid_dist_key.hpp Grid/staggered_dist_grid.hpp Grid/staggered_dist_grid_util.hpp Grid/staggered_dist_grid_copy.hpp \
         Vector/se_class3_vector.hpp  Vector/vector_dist_multiphase_functions.hpp Vector/vector_dist_comm.hpp Vector/vector_dist.hpp Vector/vector_dist_ofb.hpp Vector/Iterators/vector_dist_iterator.hpp Vector/vector_dist_key.hpp \
         config/config.h \
         example.mk \
         Decomposition/Distribution/metis_util.hpp Decomposition/Distribution/SpaceDistribution.hpp Decomposition/Distribution/parmetis_dist_util.hpp  Decomposition/Distribution/parmetis_util.hpp Decomposition/Distribution/MetisDistribution.hpp Decomposition/Distribution/ParMetisDistribution.hpp Decomposition/Distribution/DistParMetisDistribution.hpp  Decomposition/dec_optimizer.hpp SubdomainGraphNodes.hpp \
         Graph/ids.hpp Graph/dist_map_graph.hpp Graph/DistGraphFactory.hpp \
         DLB/DLB.hpp DLB/LB_Model.hpp

#testa_SOURCES = Decomposition/Domain_NN_calculator_cart_unit_test.cpp ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_vcluster/src/VCluster/VCluster.cpp
#testa_LDADD = $(LINKLIBS)
#testa_CXXFLAGS = $(OPENMP_CFLAGS) $(AM_CXXFLAGS) $(LIBHILBERT_INCLUDE) $(PETSC_INCLUDE) $(HDF5_CPPFLAGS) $(CUDA_CFLAGS) $(INCLUDES_PATH) $(PARMETIS_INCLUDE) $(METIS_INCLUDE) $(BOOST_CPPFLAGS) $(H5PART_INCLUDE) -DPARALLEL_IO  -Wno-unused-local-typedefs

lib_LIBRARIES = libofpm_pdata.a
libofpm_pdata_a_SOURCES = lib/pdata.cpp
libofpm_pdata_a_CXXFLAGS = $(OPENMP_CFLAGS) $(AM_CXXFLAGS) $(INCLUDES_PATH) $(BOOST_CPPFLAGS) -I/usr/local/include -Wno-unused-local-typedefs
libofpm_pdata_a_CFLAGS =

.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<

test: pdata
	source $(HOME)/openfpm_vars &&  cd .. && mpirun -np 3 ./src/pdata && mpirun -np 4 ./src/pdata
