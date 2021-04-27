include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

vcluster: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: vcluster

run: all
	mpirun --oversubscribe -np 3 ./vcluster

.PHONY: clean all run

clean:
	rm -f *.o *~ core vcluster

