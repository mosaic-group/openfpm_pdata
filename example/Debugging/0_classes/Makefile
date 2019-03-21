include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -g3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

se_classes: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS_SE2)

all: se_classes

run: all

.PHONY: clean all run

clean:
	rm -f *.o *~ core se_classes

