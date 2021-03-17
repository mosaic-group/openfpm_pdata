include ../../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++14 -o $@ $< $(INCLUDE_PATH)

pse_1d: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: pse_1d

run: all
	./pse_1d

.PHONY: clean all run

clean:
	rm -f *.o *~ core pse_1d

