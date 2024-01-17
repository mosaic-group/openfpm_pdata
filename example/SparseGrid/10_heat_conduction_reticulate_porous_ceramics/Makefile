include ../../example.mk
include ../../common.mk

OBJ = main.o
DATA = reactiondiffusion_imagebased_porousmedia-main.tar.gz

heat_ceram:
heat_ceram_test: OPT += -DTEST_RUN
heat_ceram_test: heat_ceram

$(DATA):
	mkdir -p output_sussman_sparse_grid_porousCeramics_1216x1016x941/
	wget https://git.mpi-cbg.de/mosaic/software/image-based-modeling/reactiondiffusion_imagebased_porousmedia/-/archive/main/reactiondiffusion_imagebased_porousmedia-main.tar.gz
	tar -xzvf reactiondiffusion_imagebased_porousmedia-main.tar.gz
	mv reactiondiffusion_imagebased_porousmedia-main/porous_catalysts/data/output_sussman_sparse_grid_porousCeramics_1216x1016x941/* output_sussman_sparse_grid_porousCeramics_1216x1016x941/
	rm -rf reactiondiffusion_imagebased_porousmedia-main/

%.o: %.cu
	$(CUDA_CC) $(CUDA_OPTIONS) -c -o $@ $< $(INCLUDE_PATH_NVCC)

heat_ceram: $(OBJ) $(DATA)
	$(CUDA_CC_LINK) -o $@ $(OBJ) $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

heat_ceram2: $(OBJ) $(DATA)
	$(CUDA_CC_LINK) -o $@ $(OBJ) $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: heat_ceram

run: heat_ceram_test
	mpirun --oversubscribe -np 1 ./heat_ceram

.PHONY: clean all run

clean:
	rm -f *.o *~ core heat_ceram

