include ../../example.mk

CC=mpic++

LDIR =

OBJ_DORD = main.o

all: md_sym

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++11 $(OPT) -o $@ $< $(INCLUDE_PATH)

md_sym: $(OBJ_DORD)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: md_sym
	source $$HOME/openfpm_vars; mpirun -np 3 ./md_sym

.PHONY: clean all run

clean:
	rm -f *.o *~ core md_sym

