CC=gcc
LD=ld
CFLAGS=-O3 -D__USE_GNU=1 -std=gnu99
LDFLAGS=-lc -ldl
SO_CFLAGS=-fpic
SO_LDFLAGS=-shared

TARGETS=malloc_fail.so

all: $(TARGETS)

malloc_fail.so: malloc_fail.o
	$(LD) $(SO_LDFLAGS) -o $@ $^

malloc_fail.o: malloc_fail.c
	$(CC) $(CFLAGS) $(SO_CFLAGS) -c $^

clean:
	rm -f *.o
	rm -f $(TARGETS)
