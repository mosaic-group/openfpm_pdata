include ../../example.mk

CC=mpic++

LDIR =
OPT=-I/Users/lschulze/gsl/include/

OBJ = main.o

cp_rounded_rectangle:
cp_rounded_rectangle_test: OPT += -DTEST_RUN
cp_rounded_rectangle_test: cp_rounded_rectangle

%.o: %.cpp ../../openfpm_numerics/src/DCPSE/MonomialBasis.hpp ../../openfpm_numerics/src/level_set/particle_cp/particle_cp.hpp
	$(CC) -O3 -g  $(OPT) -g -c --std=c++17  -o $@ $< $(INCLUDE_PATH)

cp_rounded_rectangle: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: cp_rounded_rectangle

run: cp_rounded_rectangle_test
	mpirun --oversubscribe -np 2 ./cp_rounded_rectangle

.PHONY: clean all run

clean:
	rm -f *.o *~ core cp_rounded_rectangle
