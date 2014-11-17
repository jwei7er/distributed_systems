OBJS = randomwalk.o
CC = g++
DEBUG = -g
OPTIONS = -Wall
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

all : randomwalk clean

randomwalk : $(OBJS)
	@ $(CC) $(LFLAGS) $(OBJS) -O3 -o randomwalk
	
randomwalk.o : randomwalk.cxx
	@ $(CC) $(OPTIONS) $(CFLAGS) randomwalk.cxx
	
.PHONY: clean

clean : 
	@ \rm *.o
	
tar:
	tar cfv randomwalk.tar randomwalk.cxx

