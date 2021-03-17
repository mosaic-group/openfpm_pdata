include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

OPT=

all: dem

dem_test: OPT += -DTEST_RUN
dem_test: all

%.o: %.cpp
	$(CC) -O3 $(OPT) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

dem: $(OBJ)
	$(CC) -o $@ $^ $(OPT) $(CFLAGS) $(LIBS_PATH) $(LIBS)


run: dem_test
	mpirun -np 2 ./dem

.PHONY: clean all run

clean:
	rm -f *.o *~ core dem

