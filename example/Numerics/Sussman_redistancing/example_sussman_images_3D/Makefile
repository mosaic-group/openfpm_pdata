include ../../../common.mk
include ../../../example.mk

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

example_sussman_images: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: example_sussman_images

run: all
	mpirun --oversubscribe -np 2 ./example_sussman_images

.PHONY: clean all run

clean:
	rm -f *.o *~ core example_sussman_images

