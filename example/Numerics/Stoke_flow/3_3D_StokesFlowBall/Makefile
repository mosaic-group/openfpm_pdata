include ../../../common.mk
include ../../../example.mk

OBJ = main.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

Stokes3dBall: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: Stokes3dBall

run: all
	mpirun -np 4 ./Stokes3dBall

.PHONY: clean all run

clean:
	rm -f *.o *~ core Stokes3dBall
