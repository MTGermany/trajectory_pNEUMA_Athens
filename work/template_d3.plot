#!/bin/bash

project="template_d3"
for f in `ls ${project}*gnu`; do gnuplot $f; done
