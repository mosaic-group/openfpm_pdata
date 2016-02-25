SUBDIRS := $(wildcard */.)

all clean:
	for dir in $(SUBDIRS); do \
          $(MAKE) -C $$dir $@; \
        done

clean: $(SUBDIRS)

.PHONY: all clean $(SUBDIRS)

