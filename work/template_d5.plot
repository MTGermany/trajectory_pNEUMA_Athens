#!/bin/bash

project="template_d5"
for f in `ls ${project}*gnu`; do gnuplot $f; done
