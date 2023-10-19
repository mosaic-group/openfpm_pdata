include ../../example.mk
include ../../common.mk

OBJ = main.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

Active2d: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: Active2d

run: all
	mpirun -np 4 ./Active2d 41, 1e-2, 1e-3, 1

.PHONY: clean all run

clean:
	rm -f *.o *~ core Active2d
