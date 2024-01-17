# 9_inhomogeneous_diffusion_porous_catalyst_CaCO3 and 10_heat_conduction_reticulate_porous_ceramics are excluded
# due to the large data sets included. run `make run` in each corresponding folder when required
SUBDIRS := $(shell find . -mindepth 1 -maxdepth 1 -type d ! \( -name "*10_heat_conduction_reticulate_porous_ceramics" -o -name "*9_inhomogeneous_diffusion_porous_catalyst_CaCO3" \) -name "*")

all clean run: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: all clean $(SUBDIRS)

