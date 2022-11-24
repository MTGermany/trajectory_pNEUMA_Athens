#!/bin/bash

project="template_d4"
for f in `ls ${project}*gnu`; do gnuplot $f; done
