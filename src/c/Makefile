# src/c/Makefile
CC ?= cc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -pedantic
LDFLAGS ?=

all: tri_test

tri_test: tri_tri_affine.c tri_test.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f tri_test
