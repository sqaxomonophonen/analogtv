PKGS=sdl2
CFLAGS=-Wall -O3 -ffast-math -m64 $(shell pkg-config $(PKGS) --cflags)
LINK=-lm $(shell pkg-config $(PKGS) --libs)

all: test

analogtv.o: analogtv.c analogtv.h
	$(CC) $(CFLAGS) -c analogtv.c

test.o: test.c analogtv.h
	$(CC) $(CFLAGS) -c test.c

test: test.o analogtv.o
	$(CC) $(LINK) test.o analogtv.o -o test

clean:
	rm -rf *.o test
