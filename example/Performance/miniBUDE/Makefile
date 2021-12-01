include ../../example.mk

### This is a trick to avoid "Command not found if you no not have NVCC compiler". In practice the normal C++ compiler is used
### internally the example disable with the preprocessor its code if not compiled with nvcc 
CUDA_CC=
CUDA_CC_LINK=

CC=mpic++
ifdef HIP
        CUDA_CC=hipcc
        CUDA_OPTIONS= -D__NVCC__ -D__HIP__ -DCUDART_VERSION=11000 -D__CUDACC__ -D__CUDACC_VER_MAJOR__=11 -D__CUDACC_VER_MINOR__=0 -D__CUDACC_VER_BUILD__=0
        LIBS_SELECT=$(LIBS)
        CC=hipcc
	CUDA_CC_LINK=hipcc
else
	ifdef CUDA_ON_CPU
        	CUDA_CC=mpic++ -x c++ $(INCLUDE_PATH)
        	INCLUDE_PATH_NVCC=
        	CUDA_CC_LINK=mpic++
        	CUDA_OPTIONS=-DCUDA_ON_CPU -D__NVCC__ -DCUDART_VERSION=11000
        	LIBS_SELECT=$(LIBS)
	else
        	ifeq (, $(shell which nvcc))
                	CUDA_CC=mpic++ -x c++ $(INCLUDE_PATH)
                	INCLUDE_PATH_NVCC=
                	CUDA_CC_LINK=mpic++
			LIBS_SELECT=$(LIBS)
        	else
                	CUDA_CC=nvcc -ccbin=mpic++
                	CUDA_CC_LINK=nvcc -ccbin=mpic++
			LIBS_SELECT=$(LIBS_NVCC)
        	endif
	endif
endif

CC=mpic++

LDIR =

OBJ = main.o

miniBUDE:

%.o: %.cu
	$(CUDA_CC) -Ofast -ffast-math $(CUDA_OPTIONS) $(OPT) -c --std=c++14 -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) -O3 $(OPT) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

miniBUDE: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: miniBUDE

run: miniBUDE
	mpirun --oversubscribe -np 2 ./miniBUDE

.PHONY: clean all run

clean:
	rm -f *.o *~ core miniBUDE

