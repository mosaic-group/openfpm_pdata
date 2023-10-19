include ../../example.mk
include ../../common.mk

OBJ_DORD = main.o

all: md_sym

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

md_sym: $(OBJ_DORD)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: md_sym
	mpirun --oversubscribe -np 3 ./md_sym

.PHONY: clean all run

clean:
	rm -f *.o *~ core md_sym

