include ../../example.mk
include ../../common.mk

OBJ = main.o

sync_test: OPT += -DTEST_RUN
sync_dlb_test: sync_dlb

cuda:
%.o: %.cu
	$(CUDA_CC) $(CUDA_OPTIONS) -c -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

cuda: $(OBJ)
	$(CUDA_CC_LINK) -o  $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

sync_dlb: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: cuda

run: sync_dlb_test
	mpirun --oversubscribe -np 2 ./sync_dlb

.PHONY: clean all run

clean:
	rm -f *.o *~ core sync_dlb

