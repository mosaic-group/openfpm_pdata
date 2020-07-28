include ../../example.mk

CC=mpic++

LDIR =
OPT=

SRC	= tests.cpp
OBJ = tests.o
out = sph_dlb

sph_dlb:
sph_dlb_test: OPT += -DTEST_RUN
sph_dlb_test: sph_dlb

%.o: %.cpp
	$(CC) $(OPT) -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

sph_dlb: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: sph_dlb

pretty:
	find . -regex '.*\.\(cpp\|hpp\|cc\|cxx\)' -exec clang-format -i -style=file {} \;

run: sph_dlb_test
	mpirun -np 2 ./sph_dlb

.PHONY: clean all run

clean:
	rm -f *.o *~ core sph_dlb

build:
	$(MAKE) -b pretty
	$(CC) -g -c --std=c++11 -o $(OBJ) $(SRC) -DTEST_RUN -I -Wno-deprecated-declarations   -I.  -I/home/$(USER)/opt/openfpm/openfpm_numerics/include -I/home/$(USER)/scratch/openfpm_pdata/src/config -I/home/$(USER)/scratch/openfpm_pdata/src -I/home/$(USER)/opt/openfpm/openfpm_data/include -I/home/$(USER)/opt/openfpm/openfpm_vcluster/include -I/home/$(USER)/opt/openfpm/openfpm_io/include -I/home/$(USER)/opt/openfpm/openfpm_devices/include -I/home/$(USER)/opt/openfpm/dependencies/VCDEVEL/include  -I/home/$(USER)/opt/openfpm/dependencies/METIS/include -I/home/$(USER)/opt/openfpm/dependencies/PARMETIS/include -I/home/$(USER)/opt/openfpm/dependencies/BOOST/include -I/home/$(USER)/opt/openfpm/dependencies/HDF5/include -I/home/$(USER)/opt/openfpm/dependencies/LIBHILBERT/include
	$(CC) -o $(out) $(OBJ)  -L/home/$(USER)/opt/openfpm/openfpm_devices/lib -L/home/$(USER)/opt/openfpm/openfpm_vcluster/lib -L/home/$(USER)/opt/openfpm/dependencies/VCDEVEL/lib  -L/home/$(USER)/opt/openfpm/dependencies/METIS/lib -L/home/$(USER)/opt/openfpm/dependencies/PARMETIS/lib  -L/home/$(USER)/opt/openfpm/dependencies/BOOST/lib -L/home/$(USER)/opt/openfpm/dependencies/HDF5/lib -L/home/$(USER)/opt/openfpm/dependencies/LIBHILBERT/lib    -lvcluster -lofpmmemory -lparmetis -lmetis -lboost_iostreams -lboost_program_options -lhdf5 -llibhilbert -lVc   -lrt -ldl

test:  # todo test with mpi (see run)
	./$(out) 2>&1 | tee debug/$(shell date +%FT%T | sed -e 's/:/-/g').debug

debug:
	gdb -ex run --args ./$(out)
