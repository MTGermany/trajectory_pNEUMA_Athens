#!/bin/bash

project="template_d1"
for f in `ls ${project}*gnu`; do gnuplot $f; done
