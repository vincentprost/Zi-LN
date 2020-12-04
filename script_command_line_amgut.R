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
  make_option(c("-g", "--data_generation_model"), type="character", default="norta", 
              help="data_generation_model"),
  make_option(c("-m", "--method"), type="character", default="glasso", 
              help="method"),
  make_option(c("-s", "--seed"), type="integer", default=50, 
              help="number of variables"),
  make_option(c("-o", "--output"), type="character", default="amgut.txt", 
              help="output_file"),
  make_option(c("-j", "--julia"), type="character", default="/usr/local/bin", 
              help="location of julia home")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



source("inference.R")
source("utils.R")



data(amgut1.filt)
set.seed(opt$seed)
method = opt$method
topo = opt$topology
seq_depth = "unif"

depths <- rowSums(amgut1.filt)
amgut1.filt.n  <- t(apply(amgut1.filt, 1, norm_to_total))
amgut1.filt.cs <- round(amgut1.filt.n * min(depths))

d <- ncol(amgut1.filt.cs)
n <- nrow(amgut1.filt.cs)
e <- d

set.seed(10010)
graph <- make_graph(topo, d, e)
Prec  <- graph2prec(graph)
Cor   <- cov2cor(prec2cov(Prec))

n = 500
params = get_comm_params(amgut1.filt.cs, 2, "zinegbin")
paramat <- do.call('rbind', params)
paramat <- data.frame(apply(paramat, 2, as.numeric))

X = rmvzinegbin_new(n, ks=paramat$size, munbs=paramat$munb,
   ps=paramat$pstr0, Sigma=Cor)


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
    path[[k]] = S >= threshold[k]
  }
}
if(method == "flashweave"){

  dir.create("data", showWarnings = FALSE)
  # TODO 
  julia_setup(JULIA_HOME = opt$julia)
  output_without_ext = tools::file_path_sans_ext(opt$output)

  write.csv(X, file = paste("data/", output_without_ext, "_X.csv", sep = ""))
  julia_command("using FlashWeave")
  julia_command(paste('data_path = "data/', output_without_ext, '_X.csv"', sep = ""))
  edge_names = paste0("V", 1:d)
  
  path = list()
  nlamda = 12
  threshold = lamda_path(0.7, 1e-5, nlamda)


  for(k in 1:length(threshold)){
    julia_command(paste('netw_results = learn_network(data_path, sensitive=true, heterogeneous=false, alpha = ', threshold[k] ,')', sep = ""))
    cat('save_network("data/',output_without_ext, '.edgelist", netw_results)', sep = "")
    julia_command(paste('save_network("data/',output_without_ext, '.edgelist", netw_results)', sep = ""))

    edge_list <- try(read.table(paste("data/", output_without_ext, ".edgelist", sep = "")))
    if(inherits(edge_list, "try-error")){
      S = matrix(0, d, d) 
    }
    else{
      v1 = factor(edge_list[,1], levels = edge_names)
      v2 = factor(edge_list[,2], levels = edge_names)
      S = matrix(0, d, d)
      for(i in 1:length(v1)){
        S[v1[i], v2[i]] = 1
        S[v2[i], v1[i]] = 1
      }
    }
    path[[k]] = S
  }
}


pr = precision_recall(path, graph, verbose=FALSE, plot = T)
cat(pr$AUC, method, d, topo, n, "\n", sep = "\t")
write(paste(pr$AUC, method, d, topo, n, sep = "\t"), file=opt$output, append=TRUE)


