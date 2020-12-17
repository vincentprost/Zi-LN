library(ggplot2)
library(ggh4x)


aupr.df = read.table( "output.txt", col.names = c("aupr", "method", "d", "topo", "n", "model", "sparsity"))
df_sum = data.frame(mean = numeric(), sd = numeric(), method = character(), n = character(), d = numeric(), upper_quantile = numeric(), lower_quantile = numeric(), topo = character(), topo_title = character(), dimension_title = character())


n_unique = sort(unique(aupr.df["n"])[["n"]])
d_unique = sort(unique(aupr.df["d"])[["d"]])
topo_unique  = unique(aupr.df["topo"])[["topo"]]
method_unique = unique(aupr.df["method"])[["method"]]


for(method_ in method_unique) {
  for(d_ in d_unique) {
    for(n_ in n_unique){
      for(topo_ in topo_unique) {
        print(c(method_, d_, n_ , topo_))
        sub_sb = subset(aupr.df, (method == method_) & (d == d_) & (n == n_) & (topo == topo_))
        values = sub_sb[,"aupr"]
        m = median(values)
        s = sd(values)
        #upper_quantil
        upper_quantile = quantile(values)[4]
        lower_quantile = quantile(values)[2]
        topo_title_ = "Topology"
        dimension_title_ = "p"
        #df_sum = rbind(df_sum, list(method = method_, d = d_, n = n_, topo = topo_, upper_quantile = upper_quantile, lower_quantile = lower_quantile, topo = topo, sd = s, mean = mean))
        df_sum = rbind(df_sum, list(mean = m,  sd = s, method = method_, n = toString(n_),  d = d_, upper_quantile = upper_quantile, lower_quantile = lower_quantile, 
                                    topo = topo_, topo_title = topo_title_, dimension_title = dimension_title_), stringsAsFactors=FALSE)
      }
    }
  }
}

zero_names <- list(
  "topo" = list("band" = "Band", "erdos_renyi" = "Erdos Renyi", "scale_free" = "Scale Free"),
  "d" = list("100" = "100", "300" = "300", "500" = "500"),
  "dimension_title" = list("p" = "p", "p" = "p", "p" = "p"),
  "topo_title"  = list("Topology" = "Topology", "Topology" = "Topology", "Topology" = "Topology")
)

zero_labeller <- function(variable, value){
  return(zero_names[[variable]])
}


df_sum$method <- factor(df_sum$method, levels=c("Spiec-Easi-glasso", "Spiec-Easi-mb", "sparcc", "flashweave", "glasso", "mb", "magma-glasso", "magma-mb", "ZiLN-glasso", "ZiLN-mb"))


mynamestheme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (20)), 
                      legend.title = element_text(colour = "steelblue",  face = "bold.italic", family = "Helvetica"), 
                      legend.text = element_text(face = "italic", colour="steelblue4",family = "Helvetica"), 
                      axis.title = element_text(family = "Helvetica", size = (15), colour = "steelblue4"),
                      axis.text = element_text(family = "Courier", colour = "cornflowerblue", size = (15)))



p <- ggplot(df_sum, aes(x = n, y = mean, fill = method))  + facet_nested( dimension_title + d ~ topo_title + topo, labeller = zero_labeller) + 
  geom_bar(stat="identity", position=position_dodge(0.8), width = 0.8) +
  scale_x_discrete(limits= c("50", "100", "300", "1000")) + 
  scale_fill_manual(labels = c("Spiec-Easi-gLasso", "Spiec-Easi-MB", "Sparcc", "Flashweave", "Raw-gLasso", "Raw-MB", "MAGMA-gLasso", "MAGMA-MB", "ZiLN-glasso", "ZiLN-mb"), values=c("#9ecae1", "#3182bd", "#411f16", "#ab0a20", "#a1d99b", "#31a354", "#bdbdbd", "#636363", "#fdae6b", "#e6550d")) +
  xlab("n") +
  ylab("AUPR") +
  ylim(0,1) + 
  theme(text = element_text(size=14), axis.text.x = element_text(size=12), legend.position="bottom", legend.text=element_text(size=10)) +
  geom_errorbar(aes(ymin=lower_quantile, ymax=upper_quantile), width=.2,
                position=position_dodge(0.8))


pdf("aupr_fig_2.pdf", 8, 8)
print(p)
dev.off()



