CC	= g++
INCLUDE	= ../eigen
DEBUG	= -g
CPPOBJ	= quasiOneD.o
CFLAGS	= -Wall -c $(DEBUG)
LFLAGS  = -Wall $(DEBUG)

EXEC	= p1.exe

$(EXEC) : $(CPPOBJ)
	$(CC) $(LFLAGS) $(CPPOBJS) -o $(EXEC)

quasiOneD.o : quasiOneD.cpp
	$(CC) $(CFLAGS) quasiOneD.cpp

clean:
	\rm *.o p1
