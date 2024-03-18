#!/bin/bash

set -e

#args=(500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500)
args=(2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800)

for arg in "${args[@]}"
do
        ./build/lab3-seq $arg
done

for arg in "${args[@]}"
do
	./build/lab3-omp $arg
done
