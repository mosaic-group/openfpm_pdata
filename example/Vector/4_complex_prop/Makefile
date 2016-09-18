include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o
OBJ_SER = main_ser.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

vect_cp: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

vect_ser: $(OBJ_SER)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: vect_cp vect_ser

run: all
	source $$HOME/openfpm_vars && mpirun -np 2 ./vect_cp && mpirun -np 2 ./vect_ser

.PHONY: clean all run

clean:
	rm -f *.o *~ core vect

