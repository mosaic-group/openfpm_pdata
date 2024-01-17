include ../../../common.mk
include ../../../example.mk

OBJ = Bump.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

Bump: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
all: Bump
run: all
	mpirun -np 4 ./Bump ./alpha0.csv 2 1e-3 10 1e-4 0 1 0

.PHONY: clean all run

clean: 
	rm -f *.o *~ core Bump
