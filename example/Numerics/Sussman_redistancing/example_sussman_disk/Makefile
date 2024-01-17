include ../../../common.mk
include ../../../example.mk

OBJ = main.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

example_sussman_disk: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: example_sussman_disk

run: all
	mpirun --oversubscribe -np 2 ./example_sussman_disk

.PHONY: clean all run

clean:
	rm -f *.o *~ core example_sussman_disk

