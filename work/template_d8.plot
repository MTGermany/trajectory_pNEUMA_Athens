#!/bin/bash

project="template_d8"
for f in `ls ${project}*gnu`; do gnuplot $f; done
