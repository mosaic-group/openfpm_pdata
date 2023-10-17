include ../../common.mk
include ../../example.mk

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

cell: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: cell

run: all
	mpirun --oversubscribe -np 2 ./cell

.PHONY: clean all run

clean:
	rm -f *.o *~ core cell

