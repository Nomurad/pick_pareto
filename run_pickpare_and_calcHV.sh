#! /bin/sh 

python de_duplication.py
echo $1 | ./a.out
./wfg pareto.txt 1.1 1.1 1.1 1.1
