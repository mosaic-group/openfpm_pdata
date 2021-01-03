include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

multip: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: multip

run: all
	mpirun -np 2 ./multip

.PHONY: clean all run

clean:
	rm -f *.o *~ core multip

