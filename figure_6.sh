#!/bin/bash

for seed {1..10};
for d in 100 300 500; 
do for topo in band erdos_renyi scale_free; 
do for n in 50 100 300 1000; 
do for method in Spiec-Easi-glasso Spiec-Easi-mb glasso mb magma-glasso magma-mb ZiLN-glasso ZiLN-mb sparcc flashweave; 
do Rscript script_command_line.R -d $d -t $topo -n $n -m $method -s $seed -o norta_$seed.txt -g norta; 
done; done; done; done; done

cat norta_*.txt > norta.txt
rm norta_*.txt

Rscript draw_figure_6.R
