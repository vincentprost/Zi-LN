#!/bin/bash

for seed {1..10};
do for n in 100; 
do for zeros in 0 0.1 0.5 0.7 0.9;
do for model in norta PZiLN; 
do for method in Spiec-Easi-glasso Spiec-Easi-mb glasso mb magma-glasso magma-mb ZiLN-glasso ZiLN-mb sparcc flashweave; 
do Rscript script_command_line.R -d 300 -t erdos_renyi -g $model -z $zeros -y 0.98 -n 100 -m $method -s $seed -o zeros_ziln_$seed.txt; 
done; done; done; done; done

cat zeros_ziln_*.txt > zeros_ziln.txt
rm zeros_ziln_*.txt

Rscript draw_figure_3.R
