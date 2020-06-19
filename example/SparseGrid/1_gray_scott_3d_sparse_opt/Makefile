include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -mavx -g  -O3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

gray_scott_sparse_opt: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: gray_scott_sparse_opt

run: all
	mpirun -np 4 ./gray_scott_sparse_opt

.PHONY: clean all run

clean:
	rm -f *.o *~ core gray_scott_sparse_opt

