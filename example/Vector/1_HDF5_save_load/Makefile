include ../../example.mk

CC=mpic++

LDIR =
OPT=

OBJ = main.o

all: hdf5


%.o: %.cpp
	$(CC) -O3 $(OPT) -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

hdf5: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)


run: hdf5
	mpirun -np 2 ./hdf5

.PHONY: clean all run

clean:
	rm -f *.o *~ core hdf5

