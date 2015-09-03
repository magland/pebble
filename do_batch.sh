#!/bin/bash

for i in $(seq 0 511);
do
	srun bin/pebble $i &
done
