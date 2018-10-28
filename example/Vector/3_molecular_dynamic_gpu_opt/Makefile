include ../../example.mk

CC_SCOREP=scorep --nocompiler  --cuda --mpp=mpi nvcc
CC=nvcc #${CC_SCOREP}
CC_MPI=mpic++

LDIR =

OBJ_GPU = main_gpu.o
OBJ_CPU = main_cpu.o
OBJ_CPU_BEST = main_cpu_best.o

all: md_dyn_gpu md_dyn_cpu md_dyn_cpu_best

%.o: %.cu
	$(CC) -O3 -g -c -isystem=/home/i-bird/MPI/include --std=c++11 -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC_MPI) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

md_dyn_gpu: $(OBJ_GPU)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH_NVCC) $(LIBS) -L/home/i-bird/MPI/lib  -L/usr/local/cuda/lib64 -lcudart -lmpi -L/usr/local/cuda/extras/CUPTI/lib64 -lhdf5

md_dyn_cpu: $(OBJ_CPU)
	$(CC_MPI) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

md_dyn_cpu_best: $(OBJ_CPU_BEST)
	$(CC_MPI) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: all
	mpirun -np 3 ./md_dyn_gpu && mpirun -np 3 ./md_dyn_cpu && mpirun -np 3 ./md_dyn_cpu_best;

.PHONY: clean all run

clean:
	rm -f *.o *~ core md_dyn_gpu md_dyn_cpu md_dyn_cpu_best

