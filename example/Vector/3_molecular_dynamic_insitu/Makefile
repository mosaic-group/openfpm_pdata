include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

all: md_dyn

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

md_dyn: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)


run: all
	mpirun -np 3 ./md_dyn;

.PHONY: clean all run

clean:
	rm -f *.o *~ core md_dyn

