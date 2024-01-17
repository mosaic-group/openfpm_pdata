include ../../example.mk
include ../../common.mk

OBJ = main.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

periodic: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: periodic

run: all
	mpirun --oversubscribe -np 4 ./periodic

.PHONY: clean all run

clean:
	rm -f *.o *~ core periodic

