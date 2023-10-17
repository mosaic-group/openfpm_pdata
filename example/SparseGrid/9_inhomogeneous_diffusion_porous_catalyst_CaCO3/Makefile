include ../../common.mk
include ../../example.mk

OBJ = main.o
DATA = reactiondiffusion_imagebased_porousmedia-main.tar.gz

diff_CaCO3:
diff_CaCO3_test: OPT += -DTEST_RUN
diff_CaCO3_test: diff_CaCO3

$(DATA):
	mkdir -p output_sussman_maxIter6e3_CaCO3_fluidPhase_531x531x531/
	wget https://git.mpi-cbg.de/mosaic/software/image-based-modeling/reactiondiffusion_imagebased_porousmedia/-/archive/main/reactiondiffusion_imagebased_porousmedia-main.tar.gz
	tar -xzvf reactiondiffusion_imagebased_porousmedia-main.tar.gz
	mv reactiondiffusion_imagebased_porousmedia-main/porous_catalysts/data/output_sussman_CaCO3_fluidPhase_531x531x531/* output_sussman_maxIter6e3_CaCO3_fluidPhase_531x531x531/
	rm -rf reactiondiffusion_imagebased_porousmedia-main/

%.o: %.cu
	$(CUDA_CC) -O3 $(OPT) $(CUDA_OPTIONS) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH_NVCC)

diff_CaCO3: $(OBJ) $(DATA)
	$(CUDA_CC_LINK) -o $@ $(OBJ) $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

diff_CaCO32: $(OBJ) $(DATA)
	$(CUDA_CC_LINK) -o $@ $(OBJ) $(CFLAGS) $(LIBS_PATH) $(LIBS_SELECT)

all: diff_CaCO3

run: diff_CaCO3_test
	mpirun --oversubscribe -np 1 ./diff_CaCO3

.PHONY: clean all run

clean:
	rm -f *.o *~ core diff_CaCO3

