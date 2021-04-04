include ../../example.mk

### This is a trick to avoid "Command not found if you no not have NVCC compiler". In practice the normal C++ compiler is used
### internally the example disable with the preprocessor its code if not compiled with nvcc 
CUDA_CC=
CUDA_CC_LINK=
ifdef CUDA_ON_CPU
	CUDA_CC=mpic++ -x c++ $(INCLUDE_PATH)
        INCLUDE_PATH_NVCC=
        CUDA_CC_LINK=mpic++
else
	ifeq (, $(shell which nvcc))
        	CUDA_CC=mpic++ -x c++ $(INCLUDE_PATH)
        	INCLUDE_PATH_NVCC=
       		CUDA_CC_LINK=mpic++
	else
        	CUDA_CC=nvcc -ccbin=mpic++
        	CUDA_CC_LINK=nvcc -ccbin=mpic++
	endif
endif

CC=mpic++

LDIR=
OPT=

OBJ = main.o

all: md_dyn
test: md_dyn_test
md_dyn_test: OPT += -DTEST_RUN
md_dyn_test: md_dyn

%.o: %.cu
	$(CUDA_CC) $(OPT) -O3 -g -c  --std=c++14 -o $@ $< $(INCLUDE_PATH_NVCC)

md_dyn: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(LIBS_PATH) $(LIBS)


run: md_dyn_test
	mpirun --oversubscribe -np 3 ./md_dyn;

.PHONY: clean all run

clean:
	rm -f *.o *~ core md_dyn md_dyn_expr md_dyn_vl

