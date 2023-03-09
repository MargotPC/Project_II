#!/bin/bash

dir1="results"
dir2="plots"

file1="125_dynamics_0.800.xyz"
file2="125_energy_0.800.dat"
file3="initial_conf_0.800_sc.xyz"

# Check if the directories exist and create them if they don't

mv "$file1" "$dir1"/"$file1"
mv "$file2" "$dir1"/"$file2"
mv "$file3" "$dir1"/"$file3"

