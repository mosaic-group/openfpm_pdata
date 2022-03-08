include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O0 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

Active2d: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: Active2d

run: all
	mpirun -np 4 /Active2d

.PHONY: clean all run

clean:
	rm -f *.o *~ core Active2d
