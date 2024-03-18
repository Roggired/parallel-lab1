CC=gcc
CFLAGS=-O0 -Wall
LDFLAGS=-lm

SRC=./src
BUILD=./build

all: prepare lab1-seq lab1-par-1 lab1-par-less lab1-par-exact lab1-par-more lab3-seq lab3-omp

prepare:
	mkdir -p $(BUILD)

clean:
	rm -rf $(BUILD)

lab1-seq:
	$(CC) $(CFLAGS) -o $(BUILD)/lab1-seq $(SRC)/lab1.c $(LDFLAGS)

lab1-par-1:
	$(CC) $(CFLAGS) -floop-parallelize-all -ftree-parallelize-loops=1 -o $(BUILD)/lab1-par-1 $(SRC)/lab1.c $(LDFLAGS)

lab1-par-less:
	$(CC) $(CFLAGS) -floop-parallelize-all -ftree-parallelize-loops=12 -o $(BUILD)/lab1-par-less $(SRC)/lab1.c $(LDFLAGS)

lab1-par-exact:
	$(CC) $(CFLAGS) -floop-parallelize-all -ftree-parallelize-loops=24 -o $(BUILD)/lab1-par-exact $(SRC)/lab1.c $(LDFLAGS)

lab1-par-more:
	$(CC) $(CFLAGS) -floop-parallelize-all -ftree-parallelize-loops=36 -o $(BUILD)/lab1-par-more $(SRC)/lab1.c $(LDFLAGS)

lab3-seq:
	$(CC) $(CFLAGS) -o $(BUILD)/lab3-seq $(SRC)/lab3.c $(LDFLAGS)

lab3-omp:
	$(CC) $(CFLAGS) -fopenmp -o $(BUILD)/lab3-omp $(SRC)/lab3.c $(LDFLAGS)