include ../../example.mk
include ../../common.mk

OBJ = main.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

grid: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: grid

run: all
	mpirun --oversubscribe -np 2 ./grid

.PHONY: clean all run

clean:
	rm -f *.o *~ core grid

