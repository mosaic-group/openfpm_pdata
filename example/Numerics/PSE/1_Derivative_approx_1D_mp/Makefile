# This example does not work using clang
# Eliminate the comment to activate it

include ../../../example.mk

CC=mpic++

LDIR =

OBJ_128 = main_float128.o

%.o: %.cpp
	$(CC) -O3 -c --std=c++11 -o $@ $< $(INCLUDE_PATH)

# pse_1d_128: $(OBJ_128)
# 	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS) -lquadmath

#all: pse_1d_128

run: #all
	source $$HOME/openfpm_vars;# ./pse_1d_128

#.PHONY: clean all run

clean:
	rm -f *.o *~ core pse_1d_128

