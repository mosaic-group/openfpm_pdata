SUBDIRS := $(shell find . -mindepth 1 -maxdepth 1 -type d ! -name "*Bump" -name "*")

all clean run: $(SUBDIRS)

$(SUBDIRS):
	$(MAKE) -C $@ $(MAKECMDGOALS)

.PHONY: all clean $(SUBDIRS)

