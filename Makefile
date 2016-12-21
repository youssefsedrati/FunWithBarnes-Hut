CC=gcc
CFLAGS= -std=c11 -Wall -O3 -mpclmul
LDFLAGS= -std=c11 -Wall -lm
EXEC=tests
SRC=Particle.c Cell.c Quadtree.c Morton.c
OBJ=$(SRC:.c=.o)

all: $(EXEC)

tests: tests.o $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)

.PHONY: clean mrproper

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
