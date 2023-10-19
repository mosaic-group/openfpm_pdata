include ../../../common.mk
include ../../../example.mk

OBJ = Ellip.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

Ellipsoid: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
all: Ellipsoid

run: all
	mpirun -np 4 ./Ellipsoid 31 4 1.5

.PHONY: clean all run

clean: 
	rm -f *.o *~ core Ellipsoid
