include ../../example.mk

CC=mpic++

LDIR =

OBJ = main.o

%.o: %.cpp
	$(CC) -O3 -g3 -c --std=c++11 -fext-numeric-literals  -o $@ $< $(INCLUDE_PATH)

conv_p: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS) -lquadmath

all: conv_p

run: all
	source $$HOME/openfpm_vars; ./conv_p

.PHONY: clean all run

clean:
	rm -f *.o *~ core conv_p

