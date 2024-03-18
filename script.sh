#!/bin/bash

set -e

args=(6000, 220400, 434800, 649200, 863600, 1078000, 1292400, 1506800, 1721200, 1935600, 2150000)

for arg in "${args[@]}"
do
	./build/lab3-seq $arg
done
