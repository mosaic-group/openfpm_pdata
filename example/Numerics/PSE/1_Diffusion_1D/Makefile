include ../../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -fext-numeric-literals -O3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

diff_1d: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: diff_1d

run: all
	source $$HOME/openfpm_vars; mpirun -np 4 ./diff_1d

.PHONY: clean all

clean:
	rm -f *.o *~ core diff_1d

