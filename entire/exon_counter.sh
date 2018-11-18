#!/bin/bash
n=1
while [ $n -le 21 ]
do
	cut -f10 chr$n.csv > chr$n\_exon_counts.txt
	let n=$n+1
done
