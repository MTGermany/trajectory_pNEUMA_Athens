#!/bin/bash

project="template_d6"
for f in `ls ${project}*gnu`; do gnuplot $f; done
