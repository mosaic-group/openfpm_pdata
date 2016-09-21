include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o
OBJ_SER = main_ser.o

all: vect_cp vect_ser

%.o: %.cpp
	$(CC) -O0 -g -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

vect_cp: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

vect_ser: $(OBJ_SER)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

run: all
	source $$HOME/openfpm_vars && mpirun -np 2 ./vect_cp && mpirun -np 2 ./vect_ser

.PHONY: clean all run

clean:
	rm -f *.o *~ core vect_cp vect_ser

