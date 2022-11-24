#!/bin/bash

project="20181024_d4_0900_0930"
for f in `ls ${project}*gnu`; do gnuplot $f; done
