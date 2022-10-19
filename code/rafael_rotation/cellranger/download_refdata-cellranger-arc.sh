#!/bin/bash

cwd=$(pwd)

cd $MYSCRATCH

mkdir refdata-cellranger-arc
cd refdata-cellranger-arc

curl -O https://cf.10xgenomics.com/supp/cell-arc/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz

# /fastscratch/myscratch/$USER/refdata/refdata-cellranger-arc-GRCh38-2020-A-2.0.0.tar.gz

cd $cwd
