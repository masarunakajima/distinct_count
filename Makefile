
DIR=src
CC=g++
CFLAGS=-std=c11 -Wall -g -pg -I$(DIR)/smithlab_cpp -I$(DIR)
LDFLAGS=-lsmithlab_cpp
LIBS= -L$(DIR)/smithlab_cpp -lsmithlab_cpp
objects = $(DIR)/{utils,constants}.o 
exec=count_dist

#all : target 

target:
	$(MAKE) -C $(DIR)
	$(CC) $(CFLAGS) -o $(exec) $(DIR)/dist_count.cpp $(objects) $(LIBS) 

clean:
	rm -f $(exec)
	$(MAKE) -C $(DIR) clean


