include ../../example.mk

CUDA_CC=
ifeq (, $(shell which nvcc))
        CUDA_CC=mpic++ -x c++ $(INCLUDE_PATH)
        INCLUDE_PATH_NVCC=
else
        CUDA_CC=nvcc
endif

CC_SCOREP=scorep --nocompiler  --cuda --mpp=mpi nvcc
CC=$(CUDA_CC) #${CC_SCOREP}
CC_MPI=mpic++

LDIR =
OPT=

OBJ_GPU = main_gpu.o
OBJ_CPU = main_cpu.o
OBJ_CPU_BEST = main_cpu_best.o

all: md_dyn_gpu md_dyn_cpu md_dyn_cpu_best

md_dyn_gpu_test: OPT += -DTEST_RUN
md_dyn_gpu_test: all

%.o: %.cu
	$(CC) $(OPT)  -O3 -g -c -isystem=/home/i-bird/MPI/include --std=c++11 -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC_MPI) $(OPT) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

md_dyn_gpu: $(OBJ_GPU)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS) -L/home/i-bird/MPI/lib  -L/usr/local/cuda/lib64 -lcudart -lmpi -L/usr/local/cuda/extras/CUPTI/lib64 -lhdf5

md_dyn_cpu: $(OBJ_CPU)
	$(CC_MPI) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

md_dyn_cpu_best: $(OBJ_CPU_BEST)
	$(CC_MPI) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: md_dyn_gpu_test
	mpirun -np 3 ./md_dyn_gpu && mpirun -np 3 ./md_dyn_cpu && mpirun -np 3 ./md_dyn_cpu_best;

.PHONY: clean all run

clean:
	rm -f *.o *~ core md_dyn_gpu md_dyn_cpu md_dyn_cpu_best

