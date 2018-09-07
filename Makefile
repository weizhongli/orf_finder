CC = g++
CCFLAGS = -O2 -g
#LDFLAGS = -static -o
LDFLAGS = -o
PROGS = orf_finder

all: $(PROGS)
clean:
	rm *.o $(PROGS)

# programs
orf_finder: orf_finder.o
	$(CC) $(CCFLAGS) orf_finder.o $(LDFLAGS) orf_finder

orf_finder.o: orf_finder.c++
	$(CC) $(CCFLAGS) orf_finder.c++ -c
