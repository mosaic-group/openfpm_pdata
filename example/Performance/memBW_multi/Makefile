include ../../example.mk
include ../../common.mk

OBJ = main.o

memBW:

%.o: %.cu
	$(CUDA_CC) $(CUDA_OPTIONS) -c -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

memBW: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: memBW

run: memBW
	mpirun --oversubscribe -np 2 ./memBW

.PHONY: clean all run

clean:
	rm -f *.o *~ core memBW

