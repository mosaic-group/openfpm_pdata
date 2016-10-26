include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

vect: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: vect

run: all
	source $$HOME/openfpm_vars; mpirun -np 2 ./vect

.PHONY: clean all run

clean:
	rm -f *.o *~ core vect

