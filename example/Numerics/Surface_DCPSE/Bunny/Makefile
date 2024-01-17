include ../../../common.mk
include ../../../example.mk

OBJ = Bunny.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

Bunny: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)
all: Bunny
run: all
	mpirun -np 4 ./Bunny 0.005 1

.PHONY: clean all run

clean: 
	rm -f *.o *~ core Bunny
