
LINKLIBS = $(METIS_LIB) $(PARMETIS_LIB)  $(PTHREAD_LIBS) $(OPT_LIBS) $(BOOST_LDFLAGS) $(BOOST_IOSTREAMS_LIB) $(CUDA_LIBS)

noinst_PROGRAMS = cart_dec metis_dec dom_box vector_dist
cart_dec_SOURCES = CartDecomposition_gen_vtk.cpp ../src/lib/pdata.cpp ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp
cart_dec_CXXFLAGS = $(METIS_INCLUDE) $(PARMETIS_INCLUDE) $(CUDA_CFLAGS) $(INCLUDES_PATH) $(BOOST_CPPFLAGS) -I../src -I/usr/local/libhilbert/include -Wno-unused-function -Wno-unused-local-typedefs
cart_dec_CFLAGS = $(CUDA_CFLAGS)
cart_dec_LDADD = $(LINKLIBS) -lparmetis -lmetis

metis_dec_SOURCES = Metis_gen_vtk.cpp ../src/lib/pdata.cpp ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp
metis_dec_CXXFLAGS = $(METIS_INCLUDE) $(CUDA_CFLAGS) $(INCLUDES_PATH) $(BOOST_CPPFLAGS) -I../src -I/usr/local/libhilbert/include -Wno-unused-function -Wno-unused-local-typedefs
metis_dec_CFLAGS = $(CUDA_CFLAGS)
metis_dec_LDADD = $(LINKLIBS) -lmetis

dom_box_SOURCES = domain_gen_vtk.cpp ../src/lib/pdata.cpp ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_devices/src/memory/PtrMemory.cpp ../openfpm_vcluster/src/VCluster.cpp ../openfpm_devices/src/Memleak_check.cpp
dom_box_CXXFLAGS = $(METIS_INCLUDE) $(CUDA_CFLAGS) $(INCLUDES_PATH) $(BOOST_CPPFLAGS) -I../src -I/usr/local/libhilbert/include -Wno-unused-function -Wno-unused-local-typedefs
dom_box_CFLAGS = $(CUDA_CFLAGS)
dom_box_LDADD = $(LINKLIBS)

vector_dist_SOURCES = vector.cpp ../openfpm_devices/src/memory/HeapMemory.cpp ../openfpm_vcluster/src/VCluster.cpp ../openfpm_devices/src/memory/PtrMemory.cpp
vector_dist_CXXFLAGS = $(PARMETIS_INCLUDE) $(METIS_INCLUDE) $(CUDA_CFLAGS) $(INCLUDES_PATH) $(HDF5_CPPFLAGS) $(BOOST_CPPFLAGS) -I../src -I/usr/local/libhilbert/include -Wno-unused-function -Wno-unused-local-typedefs
vector_dist_CFLAGS = $(CUDA_CFLAGS)
vector_dist_LDADD = $(LINKLIBS) -lparmetis -lmetis

.cu.o :
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<


#### Rule to make images

images : cart_dec metis_dec dom_box
	rm -rf CartDecomposition
	rm -rf Metis
	rm -rf Vector
	mkdir CartDecomposition
	mkdir Metis
	mkdir Vector
	mpirun -np 4 ./cart_dec
	mpirun -np 3 ./vector_dist
	./dom_box
	./metis_dec
	mkdir -p generated
	mkdir -p static
	pvbatch Domain_decomposition.py
	pvbatch vector_dist.py
	pvbatch mooving_particles.py
	pvbatch mooving_particles_prc.py
	pvbatch vector_scal_vect.py
	pvbatch vector_particles.py
	pvbatch particles_maps.py
	pvbatch vector_ghost.py
	dot -Tsvg openfpm.dot -o generated/openfpm.svg
	avconv -i generated/particles_mooving.ogv -f mp4 generated/particles_mooving.mp4
	avconv -i generated/particles_mooving_prc.ogv -f mp4 generated/particles_mooving_prc.mp4
	avconv -i generated/particles_mooving.ogv -f webm generated/particles_mooving.webm
	avconv -i generated/particles_mooving_prc.ogv -f webm generated/particles_mooving_prc.webm
