include ../../../common.mk
include ../../../example.mk

OBJ = CircLap.o
OBJ2=CircPoisson.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

Circ: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
CircPoisson: $(OBJ2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: Circ CircPoisson

run: all
	mpirun -np 4 ./Circ 1 2 4.9 && mpirun -np 4 ./CircPoisson 256 2 2

.PHONY: clean all run

clean: 
	rm -f *.o *~ core Circ CircPoisson
