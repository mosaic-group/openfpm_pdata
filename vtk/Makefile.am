
LINKLIBS = $(METIS_LIB) $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_IOSTREAMS_LIB) $(CUDA_LIBS)

noinst_PROGRAMS = cart_dec metis_dec dom_box
cart_dec_SOURCES = CartDecomposition_gen_vtk.cpp ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp
cart_dec_CXXFLAGS = $(CUDA_CFLAGS) $(INCLUDES_PATH) $(METIS_INCLUDE) $(BOOST_CPPFLAGS) -I../src -Wno-unused-function
cart_dec_CFLAGS = $(CUDA_CFLAGS)
cart_dec_LDADD = $(LINKLIBS) -lmetis

metis_dec_SOURCES = Metis_gen_vtk.cpp ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp
metis_dec_CXXFLAGS = $(CUDA_CFLAGS) $(INCLUDES_PATH) $(METIS_INCLUDE) $(BOOST_CPPFLAGS) -I../src -Wno-unused-function
metis_dec_CFLAGS = $(CUDA_CFLAGS)
metis_dec_LDADD = $(LINKLIBS) -lmetis

dom_box_SOURCES = domain_gen_vtk.cpp ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp
dom_box_CXXFLAGS = $(CUDA_CFLAGS) $(INCLUDES_PATH) $(METIS_INCLUDE) $(BOOST_CPPFLAGS) -I../src -Wno-unused-function
dom_box_CFLAGS = $(CUDA_CFLAGS)
dom_box_LDADD = $(LINKLIBS)

.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<
