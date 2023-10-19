include ../../example.mk
include ../../common.mk

OBJ = main.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

multip: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: multip

run: all
	mpirun --oversubscribe -np 2 ./multip

.PHONY: clean all run

clean:
	rm -f *.o *~ core multip

