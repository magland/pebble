#!/bin/bash

for i in $(seq 0 511);
do
	bin/pebble stage2 --channel=$i &
done
