#!/bin/bash

project="template_d7"
for f in `ls ${project}*gnu`; do gnuplot $f; done
