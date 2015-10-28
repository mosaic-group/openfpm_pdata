LINKLIBS = $(METIS_LIB) $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_IOSTREAMS_LIB) $(CUDA_LIBS)

noinst_PROGRAMS = pdata
pdata_SOURCES = main.cpp ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp
pdata_CXXFLAGS = $(CUDA_CFLAGS) $(INCLUDES_PATH) $(METIS_INCLUDE) $(BOOST_CPPFLAGS) -Wno-unused-function
pdata_CFLAGS = $(CUDA_CFLAGS)
pdata_LDADD = $(LINKLIBS) -lmetis
nobase_include_HEADERS = Decomposition/CartDecomposition.hpp Decomposition/common.hpp Decomposition/Decomposition.hpp  Decomposition/ie_ghost.hpp \
         Decomposition/nn_processor.hpp Decomposition/ie_loc_ghost.hpp Decomposition/ORB.hpp \
         Graph/CartesianGraphFactory.hpp \
         Grid/grid_dist_id.hpp Grid/grid_dist_id_iterator.hpp Grid/grid_dist_key.hpp \
         Vector/vector_dist.hpp Vector/vector_dist_iterator.hpp Vector/vector_dist_key.hpp \
         config/config.h \
         example.mk \
         metis_util.hpp dec_optimizer.hpp SubdomainGraphNodes.hpp

.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<
