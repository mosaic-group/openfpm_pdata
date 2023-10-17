include ../../common.mk
include ../../example.mk

OBJ = main.o

sph_dlb:
sph_dlb_test: OPT += -DTEST_RUN
sph_dlb_test: sph_dlb

%.o: %.cu
	$(CUDA_CC) -O3 $(OPT) $(CUDA_OPTIONS)  -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH_NVCC)

%.o: %.cpp
	$(CC) -O3 $(OPT) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

sph_dlb: $(OBJ)
	$(CUDA_CC_LINK) -o  $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

sph_dlb2: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: sph_dlb

run: sph_dlb_test
	mpirun --oversubscribe -np 2 ./sph_dlb

.PHONY: clean all run

clean:
	rm -f *.o *~ core sph_dlb

