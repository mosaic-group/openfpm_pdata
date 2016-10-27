include ../../example.mk

CC=mpic++

LDIR =

OBJ_DORD = main.o

all: complex_use

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++11 $(OPT) -o $@ $< $(INCLUDE_PATH)

complex_use: $(OBJ_DORD)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: all_test
	source $$HOME/openfpm_vars; mpirun -np 3 ./complex_use

.PHONY: clean all run all_test on_test

clean:
	rm -f *.o *~ core complex_use

