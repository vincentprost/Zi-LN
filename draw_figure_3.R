library(ggplot2)
library(ggh4x)


aupr.df = read.table( "zeros.txt", col.names = c("aupr", "method", "d", "topo", "n", "model", "sparsity"))
df_sum = data.frame(mean = numeric(), sd = numeric(), method = character(), model = character(), sparsity = character(), upper_quantile = numeric(), lower_quantile = numeric(), topo = character())

n_unique = sort(unique(aupr.df["n"])[["n"]])
d_unique = sort(unique(aupr.df["d"])[["d"]])
topo_unique  = unique(aupr.df["topo"])[["topo"]]
method_unique = unique(aupr.df["method"])[["method"]]
model_unique = unique(aupr.df["model"])[["model"]]
sparsity_unique = unique(aupr.df["sparsity"])[["sparsity"]]


for(method_ in method_unique) {
  for(n_ in n_unique) {
    for(model_ in model_unique){
      for(sparsity_ in sparsity_unique) {
        sub_sb = subset(aupr.df, (method == method_) & (sparsity == sparsity_) & (model == model_))
        values = sub_sb[,"aupr"]
        m = median(values)
        s = sd(values)
        upper_quantile = quantile(values)[4]
        lower_quantile = quantile(values)[2]
        df_sum = rbind(df_sum, list(mean = m,  sd = s, method = method_, n = n_, sparsity = sparsity_,  model = model_, upper_quantile = upper_quantile, lower_quantile = lower_quantile, 
                                    topo = topo_), stringsAsFactors=FALSE)
      }
    }
  }
}

df_sum$sparsity_title = "Sparsity"
df_sum$model_title = "Generative model"


zero_names <- list(
  "sparsity" = list("0" = "0", "0.1" = "0.1", "0.5" = "0.5", "0.7" = "0.7", "0.9" = "0.9"),
  "model" = list("norta" = "Norta", "PZiLN" = "PZiLN"),
  "sparsity_title" = list("Sparsity" = "Sparsity", "Sparsity" = "Sparsity", "Sparsity" = "Sparsity", "Sparsity" = "Sparsity", "Sparsity" = "Sparsity"),
  "model_title"  = list("Generative model" = "Generative model", "Generative model" = "Generative model")
)

zero_labeller <- function(variable, value){
  print(c(variable, value))
  return(zero_names[[variable]])
}


df_sum$method <- factor(df_sum$method, levels=c("Spiec-Easi-glasso", "Spiec-Easi-mb", "sparcc", "flashweave", "glasso", "mb", "magma-glasso", "magma-mb", "ZiLN-glasso", "ZiLN-mb"))


mynamestheme <- theme(plot.title = element_text(family = "Helvetica", face = "bold", size = (20)), 
                      legend.title = element_text(colour = "steelblue",  face = "bold.italic", family = "Helvetica"), 
                      legend.text = element_text(face = "italic", colour="steelblue4",family = "Helvetica"), 
                      axis.title = element_text(family = "Helvetica", size = (15), colour = "steelblue4"),
                      axis.text = element_text(family = "Courier", colour = "cornflowerblue", size = (15)))


bar_width = 1.0
p <- ggplot(df_sum, aes(x = n, y = mean, fill = method))  + facet_nested( model_title + model ~ sparsity_title + sparsity, labeller = zero_labeller) + 
  geom_bar(stat="identity", position=position_dodge(bar_width), width = bar_width) +
  scale_fill_manual(labels = c("Spiec-Easi-gLasso", "Spiec-Easi-MB", "Sparcc", "Flashweave", "Raw-gLasso", "Raw-MB", "MAGMA-gLasso", "MAGMA-MB", "ZiLN-glasso", "ZiLN-mb"), values=c("#9ecae1", "#3182bd", "#411f16", "#ab0a20", "#a1d99b", "#31a354", "#bdbdbd", "#636363", "#fdae6b", "#e6550d")) +
  xlab("p = 300, n = 100") +
  ylab("AUPR") +
  ylim(0,1) +
  scale_x_discrete(limits= c(50, 100, 150)) + 
  theme(text = element_text(size=14), legend.position="bottom", axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.text=element_text(size=10)) +
  geom_errorbar(aes(ymin=lower_quantile, ymax=upper_quantile), width=.2,
                position=position_dodge(bar_width))



pdf("aupr_fig3.pdf", 8, 6)
print(p)
dev.off()



