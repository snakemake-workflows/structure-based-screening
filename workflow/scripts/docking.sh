#!/usr/bin/env bash

recID="$4"
database="$5"
dataset="$6"
name="$7"
i="$8"
out_dir="$9"

cd ${out_dir}/${recID}/${dataset}
cp $1 .
cp $2 .
cp $3 .
srun vinalc --recList ${recID}.txt --ligList ${database}_${dataset}_${name}_${i}.txt --geoList ${recID}_grid.txt
cd -
