#!/bin/bash

for f in `ls 20181024_d7_0900_0930*gnu`; do gnuplot $f; done
