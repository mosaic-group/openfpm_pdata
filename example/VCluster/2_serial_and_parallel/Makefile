include ../../common.mk
include ../../example.mk

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -g -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

serial_parallel: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: serial_parallel

run: all
	mpirun --oversubscribe -np 3 ./serial_parallel

.PHONY: clean all run

clean:
	rm -f *.o *~ core serial_parallel

