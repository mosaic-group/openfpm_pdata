include ../../example.mk
include ../../common.mk

OBJ = main.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

ghost: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: ghost

run: all
	mpirun --oversubscribe -np 2 ./ghost

.PHONY: clean all run

clean:
	rm -f *.o *~ core ghost

