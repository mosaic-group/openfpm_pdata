
LINKLIBS = $(METIS_LIB) $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_IOSTREAMS_LIB) $(CUDA_LIBS)

bin_PROGRAMS = pdata
pdata_SOURCES = main.cpp ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster.cpp ../openfpm_data/src/Memleak_check.cpp
pdata_CXXFLAGS = $(CUDA_CFLAGS) $(INCLUDES_PATH) $(METIS_INCLUDE) $(BOOST_CPPFLAGS)
pdata_CFLAGS = $(CUDA_CFLAGS)
pdata_LDADD = $(LINKLIBS) -L/usr/lib64/nvidia-bumblebee/ -lmetis

.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<
