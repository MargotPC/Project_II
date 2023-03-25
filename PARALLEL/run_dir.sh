#!/bin/bash

dir1="results"
dir2="plots"

# Check if the directories exist and create them if they don't

if [ ! -d "$dir1" ]; then
  mkdir "$dir1"
fi

if [ ! -d "$dir2" ]; then
  mkdir "$dir2"
fi


