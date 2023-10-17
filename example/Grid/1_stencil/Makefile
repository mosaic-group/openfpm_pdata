include ../../common.mk
include ../../example.mk

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

stencil: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: stencil

run: all
	mpirun --oversubscribe -np 3 ./stencil

.PHONY: clean all run

clean:
	rm -f *.o *~ core stencil

