all:
	$(MAKE) -C build $@

clean:
	$(MAKE) -C build $@

install:
	$(MAKE) -C build $@

pdata:
	$(MAKE) -C build $@

numerics:
	$(MAKE) -C build $@

.PHONY: all clean install
