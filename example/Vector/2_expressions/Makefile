include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

expr: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: expr

run: all
	source $$HOME/openfpm_vars; mpirun -np 3 ./expr

.PHONY: clean all run

clean:
	rm -f *.o *~ core expr

