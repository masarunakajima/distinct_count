SUB=smithlab_cpp

CC=g++

CFLAGS=-std=c11 -Wall -g -pg -I$(SUB)
LDFLAGS=-lsmithlab_cpp
LIBS= -L$(SUB) -lsmithlab_cpp



objects = utils.o constants.o 

all: smithlab local


#target: smithlab local
	#$(CC) $(CFLAGS) -o count_sym_wtr count_sym_wtr.cpp $(objects) $(LIBS) 
	
smithlab: 
	$(MAKE) -C $(SUB)

local:
	$(CC) $(CFLAGS) -c constants.cpp utils.cpp 

clean:
	@-rm -f *.o *.a *~
	$(MAKE) -C $(SUB) clean


