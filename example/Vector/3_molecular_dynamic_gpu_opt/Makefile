include ../../example.mk

CC_SCOREP=scorep --nocompiler  --cuda --mpp=mpi nvcc
CC=${CC_SCOREP}

LDIR =

OBJ = main.o

all: md_dyn

%.o: %.cu
	$(CC) -O3 -g -c -isystem=/home/i-bird/MPI/include --std=c++11 -o $@ $< $(INCLUDE_PATH_NVCC)

md_dyn: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH_NVCC) $(LIBS) -L/home/i-bird/MPI/lib  -L/usr/local/cuda/lib64 -lcudart -lmpi -L/usr/local/cuda/extras/CUPTI/lib64 -lhdf5


run: all
	mpirun -np 3 ./md_dyn && mpirun -np 3 ./md_dyn_expr && mpirun -np 3 ./md_dyn_vl;

.PHONY: clean all run

clean:
	rm -f *.o *~ core md_dyn md_dyn_expr md_dyn_vl

