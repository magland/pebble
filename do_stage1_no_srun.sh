#!/bin/bash

for i in $(seq 0 99);
do
        bin/pebble stage1 --channel=$i &
done

