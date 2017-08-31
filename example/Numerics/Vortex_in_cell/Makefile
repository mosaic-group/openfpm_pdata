include ../../../example.mk

CC=mpic++

LDIR =

OBJ_VIC_PETSC = main_vic_petsc.o

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

all: vic_petsc

vic_petsc: $(OBJ_VIC_PETSC)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SE2)
	
run: all
	mpirun -np 3 ./vic_petsc
	
.PHONY: clean all

clean:
	rm -f *.o *~ core vic

