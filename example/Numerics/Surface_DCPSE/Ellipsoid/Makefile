include example.mk

CC=mpic++

LDIR =

OBJ = Ellip.o


%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

Ellipsoid: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
all: Ellipsoid

run: all
	mpirun -np 4 /Ellipsoid

.PHONY: clean all run

clean: 
	rm -f *.o *~ core Ellipsoid
