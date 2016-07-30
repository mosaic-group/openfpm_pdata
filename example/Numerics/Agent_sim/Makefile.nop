include ../../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH) -I/home/i-bird/HDF5/include

agsim: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: agsim

.PHONY: clean all

clean:
	rm -f *.o *~ core agsim

