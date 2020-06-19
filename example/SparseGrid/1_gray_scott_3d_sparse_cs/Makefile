include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

gray_scott_sparse_cs: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: gray_scott_sparse_cs

run: all
	mpirun -np 4 ./gray_scott_sparse_cs

.PHONY: clean all run

clean:
	rm -f *.o *~ core gray_scott_sparse_cs

