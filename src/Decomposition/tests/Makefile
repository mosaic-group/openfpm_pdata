CC=mpic++

LDIR =
OPT=

SRC	= main.cpp
OBJ = main.o
out = out_obj

out_obj:
out_obj_test: OPT += -DTEST_RUN
out_obj_test: out_obj
mpi_nprocs  := 3

%.o: %.cpp
	$(CC) $(OPT) -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

out_obj: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: out_obj

pretty:
	find . -regex '.*\.\(cpp\|hpp\|cc\|cxx\)' -exec clang-format -i -style=file {} \;

run: out_obj_test
	mpirun -np $(mpi_nprocs) --oversubscribe ./out_obj

.PHONY: clean all run

clean:
	rm -f *.o *~ core out_obj

build:
	$(CC) -g -c --std=c++11 -o $(OBJ) $(SRC) -DSE_CLASS1 -DTEST_RUN -I -Wno-deprecated-declarations   -I.  -I/home/$(USER)/opt/openfpm/openfpm_numerics/include -I/home/$(USER)/scratch/openfpm_pdata/src/config -I/home/$(USER)/scratch/openfpm_pdata/src -I/home/$(USER)/opt/openfpm/openfpm_data/include -I/home/$(USER)/scratch/openfpm_vcluster/src -I/home/$(USER)/opt/openfpm/openfpm_io/include -I/home/$(USER)/opt/openfpm/openfpm_devices/include -I/home/$(USER)/opt/openfpm/dependencies/VCDEVEL/include  -I/home/$(USER)/opt/openfpm/dependencies/METIS/include -I/home/$(USER)/opt/openfpm/dependencies/PARMETIS/include -I/home/$(USER)/opt/openfpm/dependencies/BOOST/include -I/home/$(USER)/opt/openfpm/dependencies/HDF5/include -I/home/$(USER)/opt/openfpm/dependencies/LIBHILBERT/include
	$(CC) -o $(out) $(OBJ)  -L/home/$(USER)/opt/openfpm/openfpm_devices/lib -L/home/$(USER)/opt/openfpm/openfpm_vcluster/lib -L/home/$(USER)/opt/openfpm/dependencies/VCDEVEL/lib  -L/home/$(USER)/opt/openfpm/dependencies/METIS/lib -L/home/$(USER)/opt/openfpm/dependencies/PARMETIS/lib  -L/home/$(USER)/opt/openfpm/dependencies/BOOST/lib -L/home/$(USER)/opt/openfpm/dependencies/HDF5/lib -L/home/$(USER)/opt/openfpm/dependencies/LIBHILBERT/lib    -lvcluster -lofpmmemory -lparmetis -lmetis -lboost_iostreams -lboost_program_options -lhdf5 -llibhilbert -lVc   -lrt -ldl

test:  # todo test with mpi (see run)
	mpirun -np $(mpi_nprocs) --oversubscribe ./$(out) 2>&1 | tee debug/$(shell date +%FT%T | sed -e 's/:/-/g').debug

debug:
	mpirun -np $(mpi_nprocs) --oversubscribe konsole -e gdb -ex run --args ./$(out)

redebug:
	$(MAKE) -B build
	$(MAKE) -B debug
