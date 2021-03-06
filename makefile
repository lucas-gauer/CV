CC = g++
CFLAGS = -std=c++11 -g -Wall -Wextra -Wshadow -O3
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

rm:
	$(RM) *.pgm
	$(RM) *.ppm

run:	$(MAIN)
	./$(MAIN)