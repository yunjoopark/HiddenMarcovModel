#!/bin/bash

n=1
while [ $n -le 22 ]
do
	grep "chr$n," exon.csv > chr$n.csv
	let n=$n+1
done
