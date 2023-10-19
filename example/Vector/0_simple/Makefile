include ../../example.mk
include ../../common.mk

OBJ = main.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

vect: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: vect

run: all
	mpirun --oversubscribe -np 2 ./vect

debug:
	../../../gdbgui/launch_mpi_debugger 2 ./vect

.PHONY: clean all run

clean:
	rm -f *.o *~ core vect

