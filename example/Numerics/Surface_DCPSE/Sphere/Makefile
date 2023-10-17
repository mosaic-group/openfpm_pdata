include ../../../common.mk
include ../../../example.mk

OBJ = Sph.o
OBJ2 = SphBench.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

SphereE: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

SphereBench: $(OBJ2)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)


all: SphereE SphereBench 

run: all
	mpirun -np 4 ./SphereE

.PHONY: clean all run

clean: 
	rm -f *.o *~ core SphereE SphereBench 