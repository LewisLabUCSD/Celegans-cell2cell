#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
source activate cell2cell
for i in {1..10}; do
  for j in {1..10}; do
    python  ${DIR}"/genetic_algorithm.py" &
  done
  wait
done 2>/dev/null
conda deactivate
