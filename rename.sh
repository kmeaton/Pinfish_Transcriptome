#!/bin/sh

# KM Eaton, Auburn University, 2021
# Code associated with Eaton et al. 2021 Frontiers in Ecology and Evolution
# This code was run on the Alabama Supercomputer.

sed -r 's/^(>[A-Za-z_0-9]+)\s.+$/\1/g' trinity_assembly_cdhit.fa > trinity_assembly_cdhit_truncatednames.fa
