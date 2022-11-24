#!/bin/bash

project="template_d2"
for f in `ls ${project}*gnu`; do gnuplot $f; done
