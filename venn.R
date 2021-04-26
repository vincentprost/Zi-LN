library(VennDiagram)
library(phyloseq)
library(SpiecEasi)
library(huge)
library(rMAGMA)
library(JuliaCall)

setwd("/Users/vincentprost/Documents/These/Plos/R_package/Zi-LN")

source("utils/magma_fix.R")
source("inference.R")


load("data/ll_deep.rda")
counts_el = counts[,apply(counts > 0, 2, sum) / dim(counts)[1] > 0.2]
taxmat_el = taxmat[apply(counts > 0, 2, sum) / dim(counts)[1] > 0.2,]
counts_el = as.matrix(counts_el)

JULIA_HOME = "/Applications/Julia-1.0.app/Contents/Resources/julia/bin/"
julia_setup(JULIA_HOME = JULIA_HOME)

target_edge_number = 1200
edge_number = 0

### ZiLN

Z = infer_Z(as.matrix(counts_el))

compute_graph_with_target_number_of_edges_mb = function(X, target_edge_number, tol = 1) {
  l_min = 0
  l = 0.5
  l_max = 1
  res <- do.call(huge::huge, c(lambda = l, list(x = X, method = "mb", verbose = T)))
  edge_number = sum(res$path[[1]] != 0) /2
  while(edge_number > target_edge_number + tol || edge_number < target_edge_number - tol){
    if(edge_number > target_edge_number + tol) {
      l_ = (l_max + l) / 2
      l_min = l
      l = l_
    }
    else{
      l_ = (l_min + l) / 2
      l_max = l
      l = l_
    }
    print(l)
    res <- do.call(huge::huge, c(lambda = l, list(x = X, method = "mb", verbose = T)))
    edge_number = sum(res$path[[1]] != 0) /2
  }
  return(res$path[[1]])
}

graph.zi = compute_graph_with_target_number_of_edges_mb(Z, 1200)

X_se = t(SpiecEasi::clr(counts_el + 1, 1))
graph.se = compute_graph_with_target_number_of_edges_mb(X_se, 1200)



########### Flashweave
flashweave_wrapper = function(X, alpha = 0.01) {
  write.csv(X, file = "/Users/vincentprost/Documents/These/Plos/script_plos/data/data_mat.csv")
  
  julia_command("using FlashWeave")
  julia_command(paste('data_path = "/Users/vincentprost/Documents/These/Plos/script_plos/data/data_mat.csv"', sep = ""))
  
  julia_command(paste('netw_results = learn_network(data_path, sensitive=true, heterogeneous=false, alpha = ',alpha , ')', sep = ""))
  julia_command(paste('save_network("network.edgelist", netw_results)', sep = ""))
  
  edge_list = read.table(paste("network.edgelist", sep = ""))
  
  d = dim(counts_el)[2]
  node_names = colnames(counts_el)
  
  v1 = factor(edge_list[,1], levels = node_names)
  v2 = factor(edge_list[,2], levels = node_names)
  edges = abs(edge_list[,3])
  S_flashweave = matrix(0, d, d)
  for(i in 1:length(v1)){
    val = edges[i]
    S_flashweave[v1[i], v2[i]] = val
    S_flashweave[v2[i], v1[i]] = val
  }
  return(S_flashweave)
}


compute_graph_with_target_number_of_edges_corr = function(S, target_edge_number, tol = 1) {
  l_min = 0
  l = 0.5
  l_max = 1
  diag(S) = 0
  edge_number = sum(S > l) /2
  while(edge_number > target_edge_number + tol || edge_number < target_edge_number - tol){
    if(edge_number > target_edge_number + tol) {
      l_ = (l_max + l) / 2
      l_min = l
      l = l_
    }
    else{
      l_ = (l_min + l) / 2
      l_max = l
      l = l_
    }
    print(c(l, edge_number))
    edge_number = sum(S > l) /2
    if(l == 0){
      return(S > l)
    }
  }
  return(S > l)
}

S_flashweave = flashweave_wrapper(counts_el, 0.01)
graph.flashweave = compute_graph_with_target_number_of_edges_corr(S_flashweave, target_edge_number)

########## Sparcc 

scc = sparcc(counts_el)
graph.sparcc = compute_graph_with_target_number_of_edges_corr(abs(scc$Cor), target_edge_number)

######## MAGMA


res_norm <- magma.norm(data=counts_el, distrib="ZINB", seq_depth="GMPR")
Z <- res_norm$z
graph.magma = compute_graph_with_target_number_of_edges_mb(Z, 1200)


##### Draw Venn

z1 = graph.se != 0
z2 = graph.magma != 0
z3 = graph.zi != 0
z4 = graph.sparcc != 0
z5 = graph.flashweave != 0

graph = list(z1, z2, z3, z4, z5)

draw_venn = function(graph, name){
  
  intersections = double(10 + 10 + 5)
  cex = double(5 + 10 + 10 + 5 + 1)
  
  n = 5
  l = 0
  for(i in 1:(n - 1)){
    for(j in (i + 1):n){
      l = l + 1
      intersections[l] = sum(graph[[i]] & graph[[j]]) / 2
    }
  }
  
  for(i in 1:(n - 2)){
    for(j in (i + 1):(n - 1)){
      for(k in (j + 1):n) {
        l = l + 1
        intersections[l] = sum(graph[[i]] & graph[[j]] & graph[[k]]) / 2
      }
    }
  }
  
  for(i in 1:(n - 3)){
    for(j in (i + 1):(n - 2)){
      for(k in (j + 1):(n - 1)) {
        for(m in (k + 1):n){
          l = l + 1
          intersections[l] = sum(graph[[i]] & graph[[j]] & graph[[k]] & graph[[m]]) / 2
        }
      }
    }
  }
  a1 = sum(graph[[1]]) / 2
  a2 = sum(graph[[2]]) / 2
  a3 = sum(graph[[3]]) / 2
  a4 = sum(graph[[4]]) / 2
  a5 = sum(graph[[5]]) / 2
  a12345 = sum(graph[[1]] & graph[[2]] & graph[[3]] & graph[[4]] & graph[[5]]) / 2
  
  pdf(file=name, height=11, width=11) 
    plot.new()
    venn.plot <- draw.quintuple.venn(
      area1 = a1, 
      area2 = a2, 
      area3 = a3, 
      area4 = a4, 
      area5 = a5, 
      n12 = intersections[1], 
      n13 = intersections[2], 
      n14 = intersections[3], 
      n15 = intersections[4], 
      n23 = intersections[5],
      n24 = intersections[6],
      n25 = intersections[7],
      n34 = intersections[8],
      n35 = intersections[9],
      n45 = intersections[10],
      n123 = intersections[11],
      n124 = intersections[12],
      n125 = intersections[13],
      n134 = intersections[14],
      n135 = intersections[15],
      n145 = intersections[16], 
      n234 = intersections[17],
      n235 = intersections[18],
      n245 = intersections[19],
      n345 = intersections[20],
      n1234 = intersections[21],
      n1235 = intersections[22],
      n1245 = intersections[23],
      n1345 = intersections[24],
      n2345 = intersections[25],
      n12345 = a12345,
      category = c("Spiec-Easi", "Magma", "ZiLN", "Sparcc", "Flashweave"),
      fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
      cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), 
      cat.cex = 1.5,
      margin = 0.07,
      ind = TRUE
    );
    dev.off()
    return(c(a1, a2, a3, a4, a5, intersections, a12345))
}


draw_venn(graph, "venn.pdf")

