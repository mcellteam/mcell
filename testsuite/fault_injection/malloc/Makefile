CC=gcc
LD=ld
CFLAGS=-O3 -D__USE_GNU=1 -std=gnu99
CFLAGS_GDB=-DMALLOC_INJECTOR_BREAK_TO_GDB=1
LDFLAGS=-lc -ldl
SO_CFLAGS=-fpic
SO_LDFLAGS=-shared

TARGETS=malloc_fail.so malloc_fail_gdb.so

all: $(TARGETS)

%.so: %.o
	$(LD) $(SO_LDFLAGS) $(LDFLAGS) -o $@ $^

malloc_fail.o: malloc_fail.c
	$(CC) $(SO_CFLAGS) $(CFLAGS) -c $^

malloc_fail_gdb.o: malloc_fail.c
	$(CC) $(SO_CFLAGS) $(CFLAGS) $(CFLAGS_GDB) -o $@ -c $^

clean:
	rm -f *.o
	rm -f $(TARGETS)
