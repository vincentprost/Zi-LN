library(SpiecEasi)
library(JuliaCall)



JULIA_HOME = "/Applications/Julia-1.0.app/Contents/Resources/julia/bin/"
source("utils.R")


julia_setup(JULIA_HOME = JULIA_HOME)
output = "output.txt"

data(amgut1.filt)
depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

set.seed(10010)
graph <- make_graph('cluster', d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))

n = 500
params = get_comm_params(amgut1.filt.cs, 2, "zinegbin")
paramat <- do.call('rbind', params)
paramat <- data.frame(apply(paramat, 2, as.numeric))



X = rmvzinegbin_new(n, ks=paramat$size, munbs=paramat$munb,
   ps=paramat$pstr0, Sigma=Cor)


X = make_synth_data("amgut1.filt", "X", n = 500)
graph = read.table(paste("X", "_adj.tsv", sep=""), sep="\t")

se <- spiec.easi(X, method='mb', lambda.min.ratio=1e-2, nlambda=15)


output_without_ext = tools::file_path_sans_ext(opt$output)
write.csv(X, file = paste("data/", output_without_ext, "_X.csv", sep = ""))
julia_command("using FlashWeave")
julia_command(paste('data_path = "data/', output_without_ext, '_X.csv"', sep = ""))
#julia_command(paste('netw_results = learn_network(data_path, sensitive=true, heterogeneous=false, alpha = ', threshold[k] ,')', sep = ""))
julia_command(paste('netw_results = learn_network(data_path, sensitive=true, heterogeneous=false)', sep = ""))
julia_command(paste('save_network("data/',output_without_ext, '.edgelist", netw_results)', sep = ""))
edge_list = read.table(paste("data/", output_without_ext, ".edgelist", sep = ""))
edge_names = edge_names = paste0("V", 1:d)
v1 = factor(edge_list[,1], levels = edge_names)
v2 = factor(edge_list[,2], levels = edge_names)
S = matrix(0, d, d)
for(i in 1:length(v1)){
  S[v1[i], v2[i]] = 1
}



dir.create("tmp", showWarnings = FALSE)
output_without_ext = tools::file_path_sans_ext(output)

write.csv(X, file = paste("tmp/", output_without_ext, "_X.csv", sep = ""))
julia_command("using FlashWeave")
julia_command(paste('data_path = "tmp/', output_without_ext, '_X.csv"', sep = ""))
node_names = paste0("V", 1:d)

path = list()
nlamda = 12
threshold = (1:nlamda) / (nlamda  / 0.6)
threshold = lamda_path(0.7, 1e-5, nlamda)

for(k in 1:length(threshold)){
  julia_command(paste('netw_results = learn_network(data_path, sensitive=true, heterogeneous=false, alpha = ', threshold[k] ,')', sep = ""))
  julia_command(paste('save_network("tmp/', output_without_ext, '.edgelist", netw_results)', sep = ""))
  
  edge_list <- try(read.table(paste("tmp/", output_without_ext, ".edgelist", sep = "")))
  if(inherits(edge_list, "try-error")){
    S = matrix(0, d, d) 
  }else{
    v1 = factor(edge_list[,1], levels = node_names)
    v2 = factor(edge_list[,2], levels = node_names)
    S = matrix(0, d, d)
    for(i in 1:length(v1)){
      if (edge_list[i,3] != 0){
        S[v1[i], v2[i]] = 1 
        S[v2[i], v1[i]] = 1
      }
    }
  }
  path[[k]] = S
}
unlink("tmp", recursive=TRUE)

pr.se = precision_recall(se$est$path, graph, verbose=T, plot = T)
pr.flashweave = precision_recall(path, graph, verbose=T, plot = T)

pr.se$AUC
pr.flashweave$AUC





