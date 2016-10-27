include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -g3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

plot: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: plot

run: all
	source $$HOME/openfpm_vars; ./plot

.PHONY: clean all run

clean:
	rm -f *.o *~ core plot

