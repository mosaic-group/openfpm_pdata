########################################################################
# This file must not be overwritten neither by cmake, neither by Clion #
# must remain exactly as it is                                         #
########################################################################

all:
	$(MAKE) -C build $@

clean:
	$(MAKE) -C build $@

install:
	$(MAKE) -C build $@
	script/install_parallel_debugger

pdata:
	$(MAKE) -C build $@

numerics:
	$(MAKE) -C build $@

.PHONY: all clean install
