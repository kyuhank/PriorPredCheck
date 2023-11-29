#!/bin/bash

set -ex

#Rscript -e "install.packages(c('rmdformats'))"

make

cp *.RData /output/
#cp *.csv /output/
#cp *.html /output/
