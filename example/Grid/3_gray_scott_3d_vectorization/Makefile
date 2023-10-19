include ../../example.mk
include ../../common.mk


OBJ = main.o update_new.o

%.o: %.f90
	mpif90 -ffree-line-length-none -fno-range-check -fno-second-underscore  -fimplicit-none  -mavx -O3 -c -g -o $@ $<

%.o: %.cpp
	$(CC) $(OPT) -mavx -c -Wno-ignored-attributes -o  $@ $< $(INCLUDE_PATH)

gray_scott: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS)

all: gray_scott

run: all
	mpirun --oversubscribe -np 4 ./gray_scott

.PHONY: clean all run

clean:
	rm -f *.o *~ core gray_scott

