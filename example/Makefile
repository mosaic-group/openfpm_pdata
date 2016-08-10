SUBDIRS := $(wildcard */.)

all: $(SUBDIRS)

clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
        done

$(SUBDIRS):
	$(MAKE) -C $@


.PHONY: all clean  $(SUBDIRS)

