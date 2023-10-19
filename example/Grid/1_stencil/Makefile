include ../../example.mk
include ../../common.mk

OBJ = main.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

stencil: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: stencil

run: all
	mpirun --oversubscribe -np 3 ./stencil

.PHONY: clean all run

clean:
	rm -f *.o *~ core stencil

