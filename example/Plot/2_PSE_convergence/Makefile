include ../../example.mk
include ../../common.mk

#### this example work only if you have quadmath activated ###
#### (only if -lquadmath work) uncomment to reactivate 

#OBJ = main.o

#%.o: %.cpp
# 	$(CC) $(OPT) -c -o $@ $< $(INCLUDE_PATH)

#conv_p: $(OBJ)
#	$(CC) -o $@ $^ $(CFLAGS) $(LIBS_PATH) $(LIBS) -lquadmath

#all: conv_p

run: #all
#	source $$HOME/openfpm_vars;# ./conv_p

#.PHONY: clean all run

clean:
	rm -f *.o *~ core conv_p

