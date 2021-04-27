include ../../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

example_sussman_sphere: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: example_sussman_sphere

run: all
	mpirun --oversubscribe -np 2 ./example_sussman_sphere

.PHONY: clean all run

clean:
	rm -f *.o *~ core example_sussman_sphere

