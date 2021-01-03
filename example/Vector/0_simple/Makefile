include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O0 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

vect: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: vect

run: all
	mpirun -np 2 ./vect

.PHONY: clean all run

clean:
	rm -f *.o *~ core vect

