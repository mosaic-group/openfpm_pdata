include ../../example.mk

### This is a trick to avoid "Command not found if you no not have NVCC compiler". In practice the normal C++ compiler is used
### internally the example disable with the preprocessor its code if not compiled with nvcc 
CUDA_CC=
CUDA_CC_LINK=
ifdef CUDA_ON_CPU
	CUDA_CC=mpic++ -x c++ $(INCLUDE_PATH)
        INCLUDE_PATH_NVCC=
        CUDA_CC_LINK=mpic++
        CUDA_OPTIONS=
		LIBS_SELECT=$(LIBS)
else
	ifeq (, $(shell which nvcc))
        	CUDA_CC=mpic++ -x c++ $(INCLUDE_PATH)
        	INCLUDE_PATH_NVCC=
        	CUDA_CC_LINK=mpic++
		CUDA_OPTIONS=
		LIBS_SELECT=$(LIBS)
	else
        	CUDA_CC=nvcc -ccbin=mpic++
        	CUDA_CC_LINK=nvcc -ccbin=mpic++
		CUDA_OPTIONS=-use_fast_math  -arch=sm_61 -lineinfo
		LIBS_SELECT=$(LIBS_NVCC)
	endif
endif

ifeq ($(PROFILE),ON)
        CUDA_CC=scorep --nocompiler  --cuda --mpp=mpi nvcc -ccbin=mpic++
        CUDA_CC_LINK=scorep --nocompiler  --cuda --mpp=mpi nvcc -ccbin=mpic++
else
	CUDA_CC:=$(CUDA_CC)
	CUDA_CC_LINK:=$(CUDA_CC_LINK)
endif

LDIR =
OPT=

INCLUDE_PATH_NVCC:=$(INCLUDE_PATH_NVCC) -I./include 

OBJ = main.o

heat_ceram:
heat_ceram_test: OPT += -DTEST_RUN
heat_ceram_test: heat_ceram

%.o: %.cu
	$(CUDA_CC) -O3 $(OPT) $(CUDA_OPTIONS) --extended-lambda --expt-relaxed-constexpr -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) -O3 $(OPT) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

heat_ceram: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

heat_ceram2: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: heat_ceram

run: heat_ceram_test
	mpirun --oversubscribe -np 1 ./heat_ceram

.PHONY: clean all run

clean:
	rm -f *.o *~ core heat_ceram

