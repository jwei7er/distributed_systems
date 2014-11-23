OBJS = randomwalk.o
CC = g++
GCC = mpic++
DEBUG = -g
OPTIONS = -Wall
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)

all : randomwalk clean

randomwalk : $(OBJS)
	@ $(GCC) $(LFLAGS) $(OBJS) -O3 -o randomwalk

randomwalk.o : randomwalk.cxx
	@ $(GCC) $(OPTIONS) $(CFLAGS) randomwalk.cxx

.PHONY: clean

clean :
	@ \rm *.o

tar:
	tar cfv randomwalk.tar randomwalk.cxx
