include ../../../common.mk
include ../../../example.mk

OBJ = Bump.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

Bump: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
all: Bump
run: all
	mpirun -np 4 ./Bump

.PHONY: clean all run

clean: 
	rm -f *.o *~ core Bump
