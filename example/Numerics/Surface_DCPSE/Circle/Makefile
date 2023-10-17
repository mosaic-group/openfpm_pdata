include ../../../common.mk
include ../../../example.mk

OBJ = CircLap.o
OBJ2=CircPoisson.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

Circ: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
CircPoisson: $(OBJ2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: Circ CircPoisson

run: all
	mpirun -np 4 ./CircNoProj

.PHONY: clean all run

clean: 
	rm -f *.o *~ core Circ CircPoisson
