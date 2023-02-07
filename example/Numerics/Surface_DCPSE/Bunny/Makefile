include example.mk

CC=mpic++

LDIR =

OBJ = Bunny.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

Bunny: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
all: Bunny
run: all
	mpirun -np 4 /Bunny

.PHONY: clean all run

clean: 
	rm -f *.o *~ core Bunny
