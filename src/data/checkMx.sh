#!/bin/bash
set -eu
IFS=$'\n'
cpp=$(basename $1)
py=$(basename $2)
dir=$(dirname $1)

cut -f1 $1 | sort -h > ${dir}/barcodes.$cpp
cut -f1 $2 | sort -h > ${dir}/barcodes.$py

diff=$(diff -q ${dir}/barcodes.$cpp ${dir}/barcodes.$py)
if [[ "$diff" != "" ]]; then
	echo "Barcodes are not identical!"
	exit 1
fi

while read bx; do
	grep "${bx}" $1 | cut -f2 | tr " " "\n" | sort -n > ${dir}/${bx}.${cpp}
	grep "${bx}" $2 | cut -f2 | tr " " "\n" | sort -n > ${dir}/${bx}.${py}
	diff=$(diff -q <(grep . ${dir}/${bx}.${cpp} ) <(grep . ${dir}/${bx}.${py} ))
	if [[ "$diff" != "" ]]; then
		echo "${bx} has different minimizers!"
		exit 1
	else
		rm ${dir}/${bx}.${cpp} ${dir}/${bx}.${py}
	fi
done < ${dir}/barcodes.$cpp
rm ${dir}/barcodes.$cpp ${dir}/barcodes.$py
