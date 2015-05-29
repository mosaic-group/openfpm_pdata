
LINKLIBS = $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_IOSTREAMS_LIB) $(CUDA_LIBS)

bin_PROGRAMS = pdata
pdata_SOURCES = main.cpp ../../OpenFPM_devices/src/memory/HeapMemory.cpp ../../OpenFPM_devices/src/memory/PtrMemory.cpp ../../OpenFPM_vcluster/src/VCluster.cpp ../../OpenFPM_data/src/Memleak_check.cpp
pdata_CXXFLAGS = $(CUDA_CFLAGS) $(INCLUDES_PATH)
pdata_CFLAGS = $(CUDA_CFLAGS)
pdata_LDADD = $(LINKLIBS) -L/usr/lib64/nvidia-bumblebee/ -L../../metis_install/lib  -lmetis

.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<
