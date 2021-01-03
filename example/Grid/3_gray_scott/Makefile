include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

gray_scott: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: gray_scott

run: all
	mpirun -np 4 ./gray_scott

.PHONY: clean all run

clean:
	rm -f *.o *~ core gray_scott

