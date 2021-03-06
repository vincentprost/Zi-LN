#!/usr/bin/env Rscript

library(SpiecEasi)
library(stats)
library(rMAGMA)
library(QUIC)
library(huge)
library(JuliaCall)
library(optparse)
 
option_list = list(
  make_option(c("-d", "--dimension"), type="integer", default=50, 
              help="number of variables"),
  make_option(c("-t", "--topology"), type="character", default="band", 
            help="topology of network"),
  make_option(c("-n", "--number_of_data_points"), type="integer", default=50, 
              help="dataset file name"),
  make_option(c("-z", "--sparsity"), type="numeric", default=0.45, 
              help="proportion of zeros"),
  make_option(c("-y", "--sparsity_max"), type="numeric", default=0.9, 
              help="max of sparsity for individual OTUs"),
  make_option(c("-g", "--data_generation_model"), type="character", default="PZiLN", 
              help="data_generation_model"),
  make_option(c("-m", "--method"), type="character", default="glasso", 
              help="method"),
  make_option(c("-s", "--seed"), type="integer", default=50, 
              help="number of variables"),
  make_option(c("-o", "--output"), type="character", default="output.txt", 
              help="output_file"),
  make_option(c("-j", "--julia"), type="character", default="/usr/local/bin", 
              help="location of julia home"),
  make_option(c("-f", "--fit"), type="character", default=NULL, 
              help="data for fitting parameters")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



source("inference.R")
source("utils/utils.R")

##################################
######## Data generation #########

d = opt$dimension
topo = opt$topology
n = opt$number_of_data_points
method = opt$method
set.seed(opt$seed)

e <- d
graph <- SpiecEasi::make_graph(topo, d, e)
Prec  <- SpiecEasi::graph2prec(graph, targetCondition = 100)
Cov = SpiecEasi::prec2cov(Prec)


sparsity = opt$sparsity
pstr_max = opt$sparsity_max
t = pstr_max / sparsity - 1
pstr = pstr_max * ((1:d - 0.1) / d)^t



if(opt$data_generation_model == "PZiLN"){
  seq_depth = "TS"
  mu = runif(d, 0, 3)    
  log_Y = SpiecEasi::rmvnorm(n, mu, Sigma=Cov)
  for(i in 1:d) {
    s = sort(log_Y[,i])
    cut_index = floor(n * pstr[i]) + 1
    log_Y[log_Y[,i] < s[cut_index], i] = - Inf
  }
  Y = exp(log_Y)

  P = Y  /  apply(Y, 1, sum)
  X = matrix(0, n , d)
  sequencing_depth = qnbinom(t(runif(n)), mu=150000, size=2)

  for(i in 1:n){
    X[i,] = rmultinom(1, size =  sequencing_depth[i], P[i,])
  }

}
if(opt$data_generation_model == "norta"){
  seq_depth = "unif"
  munbs = exp(runif(d, 0, 3))
  ks = 10
  if(!is.null(opt$fit)){
    if (opt$fit == "amgut1.filt") {
      data(amgut1.filt)
      source_df = amgut1.filt
    }
    else{
      source_df = read.table(source_data, sep="\t")
    }
    depths <- rowSums(source_df)
    source_filt_n <- t(apply(source_df, 1, norm_to_total))
    source_filt_cs <- round(source_filt_n * min(depths))
    d <- ncol(source_filt_cs)
    if (n <= 0) {
      n = nrow(source_filt_cs)
    }
    e <- d
    graph <- SpiecEasi::make_graph(topo, d, e)
    Prec <- SpiecEasi::graph2prec(graph, targetCondition = 100)
    Cov <- SpiecEasi::prec2cov(Prec)
    paramat = fit_parameters(source_filt_cs)
    ks = paramat$size
    munbs = munbs=paramat$munb
    pstr = paramat$pstr0
  }
  X = rmvzinegbin_new(n, ks = ks, munbs=munbs,
                      ps= pstr, Sigma=Cov)

  while(sum(apply(X == 0, 2, sum) == n) > 0){
    print("columns full of zeros, re-sample")
    X = rmvzinegbin_new(n, ks = 10, munbs=munbs,
                      ps= pstr, Sigma=Cov)
  }
}



##################################
######## Network Inference #######

paths = list()


if(method == "Spiec-Easi-glasso") {
  se <- spiec.easi(X, method="glasso", lambda.min.ratio=0.1, nlambda=12, pulsar.select = F)
  path = se$est$path
}
if(method == "Spiec-Easi-mb") {
  se <- spiec.easi(X, method="mb", lambda.min.ratio=0.1, nlambda=12, pulsar.select = F)
  path = se$est$path
}
if(method == "glasso") {
  rho_mat = matrix(1, d, d)
  diag(rho_mat) = 0
  S = cor(X)
  lamda_max = max_off_diagonal_value(S)
  pen = lamda_path(lamda_max = lamda_max)
  a = QUIC(S, rho = rho_mat, path = pen, tol = 1e-7, msg = 0, maxIter = 100000)
  quic_path = list()
  for(k in 1:length(pen)) {
    quic_path[[k]] = a$X[,,k]
  }
  path = quic_path
}
if(method == "mb") {
  rho_mat = matrix(1, d, d)
  diag(rho_mat) = 0
  S = cor(X)
  lamda_max = max_off_diagonal_value(S)
  pen = lamda_path(lamda_max = lamda_max)
  mb = huge.mb(X, lambda = pen, verbose = F)
  path = mb$path
}
if(method == "magma-glasso") {
  magma_Stool <- magma(data = X, distrib = "ZINB", method = "glasso",  seq_depth = seq_depth, magma.select = FALSE)
  path = magma_Stool$path
}
if(method == "magma-mb") {
  magma_Stool <- magma(data = X, distrib = "ZINB", method = "mb",  seq_depth = seq_depth, magma.select = FALSE)
  path = magma_Stool$path
}
if(method == "ZiLN-glasso"){
  Z = infer_Z(X, seq_depth = seq_depth)
  rho_mat = matrix(1, d,d)
  diag(rho_mat) = 0           
  S = cor(Z)
  lamda_max = max_off_diagonal_value(S)
  pen = lamda_path(lamda_max = lamda_max)
  a = QUIC(S, rho = rho_mat, path = pen, tol = 1e-7, msg = 0, maxIter = 100000)
  quic_path_Z2 = list()
  for(k in 1:length(pen)) {
    quic_path_Z2[[k]] = a$X[,,k]
  }
  path = quic_path_Z2
}
if(method == "ZiLN-mb"){
  Z = infer_Z(X, seq_depth = seq_depth)
  S = cor(Z)
  lamda_max = max_off_diagonal_value(S)
  pen = lamda_path(lamda_max = lamda_max)
  mb_Z2 = huge.mb(Z, lambda = pen, verbose = F)
  path = mb_Z2$path
}
if(method == "sparcc"){
  scc = sparcc(X)
  S = scc$Cor
  lamda_max = max_off_diagonal_value(S)
  nlamda = 20
  threshold = (1:nlamda) / (nlamda / lamda_max)
  path = list()
  for(k in 1:length(threshold)){
    path[[k]] = abs(S) >= threshold[k]
  }
}
if(method == "flashweave"){

  dir.create("tmp", showWarnings = FALSE)
  julia_setup(JULIA_HOME = opt$julia)
  output_without_ext =basename(tools::file_path_sans_ext(opt$output))

  write.csv(X, file = paste("tmp/", output_without_ext, "_X.csv", sep = ""))
  julia_command("using FlashWeave")
  julia_command(paste('data_path = "tmp/', output_without_ext, '_X.csv"', sep = ""))
  edge_names = paste0("V", 1:d)
  
  path = list()
  nlamda = 20
  

  julia_command(paste('netw_results = learn_network(data_path, sensitive=true, heterogeneous=false)', sep = ""))
  cat('save_network("tmp/', output_without_ext, '.edgelist", netw_results)', sep = "")
  julia_command(paste('save_network("tmp/', output_without_ext, '.edgelist", netw_results)', sep = ""))

  edge_list <- read.table(paste("tmp/", output_without_ext, ".edgelist", sep = ""))
  edges = abs(edge_list[,3])
  threshold = seq(min(edges),  max(edges), (max(edges) - min(edges)) / nlamda)

  v1 = factor(edge_list[,1], levels = edge_names)
  v2 = factor(edge_list[,2], levels = edge_names)
  S = matrix(0, d, d)
  for(i in 1:length(v1)){
    w = edges[i]
    if(w != 0) {
      S[v1[i], v2[i]] = w
      S[v2[i], v1[i]] = w
    }
  }

  for(k in 1:length(threshold)){
    path[[k]] = S > threshold[k]
  }
}


pr = precision_recall(path, graph, verbose=FALSE, plot = T)
cat(pr$AUC, method, d, topo, n, "\n", sep = "\t")
write(paste(pr$AUC, method, d, topo, n, opt$data_generation_model, opt$sparsity, sep = "\t"), file=opt$output, append=TRUE)


