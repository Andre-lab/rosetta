#!/bin/bash

echo 'Building documentation for each src/doxygen.*'

mkdir html
cp doc/index.html html/

rm -rf html/core+protocols html/all_else

for file in src/Doxyfile.* ; do
    echo "Building doc: $file"
    doxygen $file
done
