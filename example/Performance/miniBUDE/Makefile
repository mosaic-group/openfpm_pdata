include ../../example.mk
include ../../common.mk

OBJ = main.o

miniBUDE:

%.o: %.cu
	$(CUDA_CC) $(CUDA_OPTIONS) -c -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

miniBUDE: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: miniBUDE

run: miniBUDE
	mpirun --oversubscribe -np 1 ./miniBUDE

.PHONY: clean all run

clean:
	rm -f *.o *~ core miniBUDE

