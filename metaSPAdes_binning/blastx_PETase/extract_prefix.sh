#!/bin/bash

files=$(find ./ -type f -name "SRR*")

for file in $files; do
  filename=$(basename "$file")
  prefix="${filename%.*}"
  echo "$prefix" >> list.txt
done
