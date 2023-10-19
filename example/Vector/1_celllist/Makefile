include ../../example.mk
include ../../common.mk

OBJ = main.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

cell: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: cell

run: all
	mpirun --oversubscribe -np 2 ./cell

.PHONY: clean all run

clean:
	rm -f *.o *~ core cell

