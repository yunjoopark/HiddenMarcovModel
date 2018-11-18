#!/bin/bash
n=1
while [ $n -le 22 ]
do
	echo `wc -l "chr$n.csv"`
	let n=$n+1
done
