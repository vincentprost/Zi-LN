for seed {1..10}; 
do for d in 100 300 500; 
do for topo in band erdos_renyi scale_free; 
do for n in 50 100 300 1000; 
do for method in Spiec-Easi-glasso Spiec-Easi-mb glasso mb magma-glasso magma-mb ZiLN-glasso ZiLN-mb sparcc flashweave; 
do Rscript script_command_line.R -d $d -t $topo -n $n -m $method -s $seed -o output_$seed.txt; 
done; done; done; done; done;

cat output_*.txt > output.txt
rm output_*.txt

Rscript draw_figure_2.R