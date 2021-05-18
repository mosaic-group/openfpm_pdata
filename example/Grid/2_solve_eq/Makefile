include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

periodic: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: periodic

run: all
	mpirun --oversubscribe -np 4 ./periodic

.PHONY: clean all run

clean:
	rm -f *.o *~ core periodic

