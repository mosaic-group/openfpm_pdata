include ../../example.mk
include ../../common.mk

OBJ = main.o

%.o: %.cpp
	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

expr: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: expr

run: all
	mpirun --oversubscribe -np 3 ./expr

.PHONY: clean all run

clean:
	rm -f *.o *~ core expr

