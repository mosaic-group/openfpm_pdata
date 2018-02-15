include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

ps_cma_es: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: ps_cma_es

run: all
	mpirun -np 2 ./ps_cma_es

.PHONY: clean all run

clean:
	rm -f *.o *~ core ps_cma_es

