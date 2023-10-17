# Surface_DCPSE and 2D_ActiveFluid are excluded due to missing command line arguments
# OdeInt is excluded due to a pending issue with OdeInt
SUBDIRS := $(shell find . -mindepth 1 -maxdepth 1 -type d ! \( -name "*Surface_DCPSE" -o -name "*2D_ActiveFluid" -o -name "*OdeInt" \) -name "*")

all clean run: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: all clean $(SUBDIRS)

