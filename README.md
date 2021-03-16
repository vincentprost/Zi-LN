# Zi-LN

This repository is the code base used for generating the experimental results of the article "A zero inflated log-normal model for inference of sparse microbial association networks". The inference method exposed in the article is implemented in the scripts [inference.R](inference.R) and [script_command_line.R](script_command_line.R) along with other methods used for benchmark.



### R package dependencies

* SpiecEasi (1.0.7)
* stats
* rMAGMA (0.93)
* QUIC (1.1)
* huge (1.3.4)
* JuliaCall (0.17.1)
* optparse (1.6.6)
* ggplot2 (3.3.2)
* ggh4x

### Julia package dependencie

* FlashWeave

### Figure 2

```bash
./draw_figure_2.sh
```

[aupr_fig_2.pdf](aupr_fig_2.pdf)

### Figure 3


```bash
./draw_figure_3.sh
```

[aupr_fig_3.pdf](https://github.com/vincentprost/Zi-LN/blob/master/aupr_fig_3.pdf)

### Figure 6

```bash
./draw_figure_6.sh
```

[aupr_fig_6.pdf](https://github.com/vincentprost/Zi-LN/blob/master/aupr_fig_6.pdf)

