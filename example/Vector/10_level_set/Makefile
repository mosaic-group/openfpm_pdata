include ../../example.mk

CC=mpic++

LDIR =
OPT=

OBJ = main.o

particle_level_set:
particle_level_set_test: OPT += -DTEST_RUN
particle_level_set_test: particle_level_set

%.o: %.cpp
	$(CC) -O3 -g  $(OPT) -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

particle_level_set: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: particle_level_set

run: particle_level_set_test
	mpirun --oversubscribe -np 2 ./particle_level_set

.PHONY: clean all run

clean:
	rm -f *.o *~ core particle_level_set

