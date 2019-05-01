CC = g++
CFLAGS = -std=c++11 -Wall -g
SRCS = main.cpp image.cpp channel.cpp
OBJS = $(SRCS:.cpp=.o)
MAIN = a.out

all:    $(MAIN)
	@echo  Image Handler has been compiled
	
$(MAIN): $(OBJS) 
	$(CC) $(CFLAGS) -o $(MAIN) $(OBJS)
	
.cpp.o:
	$(CC) $(CFLAGS) -c $<  -o $@
	
clean:
	$(RM) *.o *~ $(MAIN)
	$(RM) *.pgm
	$(RM) *.ppm
