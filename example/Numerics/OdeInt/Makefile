include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o
OBJ2 = main2.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

example_odeint: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

example_odeint2: $(OBJ2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: example_odeint example_odeint2

run: all
	mpirun -np 2 ./example_odeint ./example_odeint

.PHONY: clean all run

clean:
	rm -f *.o *~ core example_odeint example_odeint2

