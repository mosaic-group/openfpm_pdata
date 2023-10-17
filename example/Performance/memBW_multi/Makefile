include ../../common.mk
include ../../example.mk

OBJ = main.o

memBW:

%.o: %.cu
	$(CUDA_CC) -O3 -g  $(CUDA_OPTIONS) $(OPT) -c --std=c++14  -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) -g -O3 $(OPT) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

memBW: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: memBW

run: memBW
	mpirun --oversubscribe -np 2 ./memBW

.PHONY: clean all run

clean:
	rm -f *.o *~ core memBW

