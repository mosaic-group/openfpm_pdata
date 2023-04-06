include ../../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

all: sphere_redistance

sphere_redistance: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: all
	mpirun --oversubscribe -np 3 ./sphere_redistance 
	
.PHONY: clean all

clean:
	rm -f *.o *~ core sphere_redistance

