include ../../example.mk

CC=mpic++

LDIR =
OPT=-I/Users/lschulze/gsl/include/ -I/Users/lschulze/minter/src

OBJ = main.o

cp_ellipsoid:
cp_ellipsoid_test: OPT += -DTEST_RUN
cp_ellipsoid_test: cp_ellipsoid

%.o: %.cpp
	$(CC) -O3 -g  $(OPT) -g -c --std=c++17  -o $@ $< $(INCLUDE_PATH)

cp_ellipsoid: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS) -L/Users/lschulze/minter/build/ -lminter 

all: cp_ellipsoid

run: cp_ellipsoid_test
	mpirun --oversubscribe -np 2 ./cp_ellipsoid

.PHONY: clean all run

clean:
	rm -f *.o *~ core cp_ellipsoid
