#!/bin/bash

for i in $(seq 0 511);
do
	srun bin/pebble stage1 --channel=$i &
done
