include ../../example.mk
include ../../common.mk

OBJ = main.o

all: md_dyn
test: md_dyn_test
md_dyn_test: OPT += -DTEST_RUN
md_dyn_test: md_dyn

%.o: %.cu
	$(CUDA_CC) $(CUDA_OPTIONS) -c -o $@ $< $(INCLUDE_PATH_NVCC)

md_dyn: $(OBJ)
	$(CUDA_CC_LINK) -o $@ $^ $(LIBS_PATH) $(LIBS_SELECT)


run: md_dyn_test
	mpirun --oversubscribe -np 3 ./md_dyn;

.PHONY: clean all run

clean:
	rm -f *.o *~ core md_dyn md_dyn_expr md_dyn_vl

