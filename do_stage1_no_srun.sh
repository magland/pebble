#!/bin/bash

for i in $(seq 0 511);
do
        bin/pebble stage1 --channel=$i &
done

