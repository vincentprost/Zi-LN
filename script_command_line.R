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
  make_option(c("-m", "--method"), type="character", default="glasso", 
              help="method"),
  make_option(c("-s", "--seed"), type="integer", default=50, 
              help="number of variables"),
  make_option(c("-o", "--output"), type="character", default="output.txt", 
              help="output_file"),
  make_option(c("-j", "--julia"), type="character", default="usr/local/bin", 
              help="location of julia home")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)



source("inference.R")
source("utils.R")



d = opt$dimension
topo = opt$topology
n = opt$number_of_data_points
method = opt$method
set.seed(opt$seed)

e <- d
graph <- SpiecEasi::make_graph(topo, d, e)
Prec  <- SpiecEasi::graph2prec(graph, targetCondition = 100)
Cov = SpiecEasi::prec2cov(Prec)


mu = runif(d, 0, 3)        

log_Y = SpiecEasi::rmvnorm(n, mu, Sigma=Cov)
log_Y_old = log_Y

pstr = runif(d, 0, 0.9)
for(i in 1:d) {
  s = sort(log_Y[,i])
  cut_index = floor(n * pstr[i]) + 1
  log_Y[log_Y[,i] < s[cut_index], i] = - Inf
}

Y = exp(log_Y)

P = Y  /  apply(Y, 1, sum)
X = matrix(0, n , d)
s_mu = 1500000
size = 5
s_var = s_mu + s_mu^2 / size

sequencing_depth = rmvnegbin(n, matrix(s_mu * d, 1, 1), matrix(s_var, 1, 1))

for(i in 1:n){
  X[i,] = rmultinom(1, size =  sequencing_depth[i], P[i,])
}



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
  magma_Stool <- magma(data = X, distrib = "ZINB", method = "glasso",  seq_depth = "TS", magma.select = FALSE)
  path = magma_Stool$path
}
if(method == "magma-mb") {
  magma_Stool <- magma(data = X, distrib = "ZINB", method = "mb",  seq_depth = "TS", magma.select = FALSE)
  path = magma_Stool$path
}
if(method == "ZiLN-glasso"){
  Z = infer_Z(X)
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
  Z = infer_Z(X)
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
  # TODO 
  julia_setup(JULIA_HOME = opt$julia)
  write.csv(X, file = "data/X.csv")
  julia_command("using FlashWeave")
  julia_command('data_path = "data/X.csv"')
  edge_names = paste0("V", 1:d)
  
  path = list()
  nlamda = 20
  threshold = (1:nlamda) / (nlamda  / 0.95)

  for(k in 1:length(threshold)){
    julia_command(paste('netw_results = learn_network(data_path, sensitive=true, heterogeneous=false, alpha = ', threshold[k] ,')'))
    julia_command('save_network("data/network_output.gml", netw_results)')
    julia_command('save_network("data/network_output.edgelist", netw_results)')
    edge_list = read.table("data/network_output.edgelist", stringsAsFactors = F)
    v1 = factor(edge_list[,1], levels = edge_names)
    v2 = factor(edge_list[,2], levels = edge_names)
    S = matrix(0, d, d)
    for(i in 1:length(v1)){
      S[v1[i], v2[i]] = 1
    }
    path[[k]] = S
  }
}


pr = precision_recall(path, graph, verbose=FALSE)
cat(pr$AUC, method, d, topo, n, "\n", sep = "\t")
write(paste(pr$AUC, method, d, topo, n, sep = "\t"),file=opt$output,append=TRUE)





