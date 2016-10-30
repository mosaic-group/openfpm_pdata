include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

ghost: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: ghost

run: all
	source $$HOME/openfpm_vars; mpirun -np 2 ./ghost

.PHONY: clean all run

clean:
	rm -f *.o *~ core ghost

