LINKLIBS = $(LIBHILBERT_LIB)  $(METIS_LIB) $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_IOSTREAMS_LIB) $(CUDA_LIBS) $(PETSC_LIB) $(HDF5_LDFLAGS)  $(HDF5_LIBS) $(PARMETIS_LIB) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) $(BOOST_CHRONO_LIB) $(BOOST_TIMER_LIB) $(BOOST_SYSTEM_LIB)

noinst_PROGRAMS = pdata
pdata_SOURCES = main.cpp Grid/grid_dist_id_unit_test.cpp  lib/pdata.cpp test_multiple_o.cpp ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp
pdata_CXXFLAGS = $(AM_CXXFLAGS) $(LIBHILBERT_INCLUDE) $(PETSC_INCLUDE) $(HDF5_CPPFLAGS) $(CUDA_CFLAGS) $(INCLUDES_PATH) $(PARMETIS_INCLUDE) $(METIS_INCLUDE) $(BOOST_CPPFLAGS) $(H5PART_INCLUDE) -DPARALLEL_IO  -Wno-unused-local-typedefs
pdata_CFLAGS = $(CUDA_CFLAGS)
pdata_LDADD = $(LINKLIBS) -lparmetis -lmetis
nobase_include_HEADERS = Decomposition/CartDecomposition.hpp Decomposition/CartDecomposition_ext.hpp  Decomposition/common.hpp Decomposition/Decomposition.hpp  Decomposition/ie_ghost.hpp \
         Decomposition/nn_processor.hpp Decomposition/ie_loc_ghost.hpp Decomposition/ORB.hpp \
         Graph/CartesianGraphFactory.hpp \
         Grid/grid_dist_id.hpp Grid/grid_dist_id_iterator_dec.hpp Grid/grid_dist_util.hpp  Grid/grid_dist_id_iterator_sub.hpp Grid/grid_dist_id_iterator.hpp Grid/grid_dist_key.hpp Grid/staggered_dist_grid.hpp Grid/staggered_dist_grid_util.hpp Grid/staggered_dist_grid_copy.hpp \
         Vector/vector_dist_comm.hpp Vector/vector_dist.hpp Vector/vector_dist_ofb.hpp Vector/vector_dist_iterator.hpp Vector/vector_dist_key.hpp \
         config/config.h \
         example.mk \
         Decomposition/Distribution/metis_util.hpp Decomposition/Distribution/parmetis_dist_util.hpp  Decomposition/Distribution/parmetis_util.hpp Decomposition/Distribution/MetisDistribution.hpp Decomposition/Distribution/ParMetisDistribution.hpp Decomposition/Distribution/DistParMetisDistribution.hpp  dec_optimizer.hpp SubdomainGraphNodes.hpp \
         Graph/ids.hpp Graph/dist_map_graph.hpp Graph/DistGraphFactory.hpp \
         DLB/DLB.hpp

lib_LIBRARIES = libofpm_pdata.a
libofpm_pdata_a_SOURCES = lib/pdata.cpp
libofpm_pdata_a_CXXFLAGS = $(INCLUDES_PATH) $(BOOST_CPPFLAGS) -I/usr/local/include -Wno-unused-local-typedefs
libofpm_pdata_a_CFLAGS =

.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<

test: pdata
	source $(HOME)/openfpm_vars &&  cd .. && mpirun -np 3 ./src/pdata && mpirun -np 4 ./src/pdata
