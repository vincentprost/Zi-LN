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

compute_graph_with_target_number_of_edges_mb = function(X, target_edge_number, tol = 50) {
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


compute_graph_with_target_number_of_edges_corr = function(S, target_edge_number, tol = 50) {
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


#####################


graph_shuffled = list()

for(k in 1:5){
  g = graph[[k]]
  edges = as.vector(g[lower.tri(g)])
  edges.shuffle = edges[sample.int(length(edges))]
  g[lower.tri(g)] = edges.shuffle
  g[upper.tri(g)] = t(g)[upper.tri(g)]
  graph_shuffled[[k]] = g
}

draw_venn(graph_shuffled, "venn_shuffle_edges.pdf")




set.seed(0)
N = 10
areas = matrix(0, N, 31)

for(k in 1:N){
  counts_el_shuffled = counts_el
  
  for(i in 1:d){
    counts_el_shuffled[, i] = sample(counts_el_shuffled[, i])
  }
  
  Z = infer_Z(counts_el_shuffled)
  graph.zi = compute_graph_with_target_number_of_edges_mb(Z, 1200)
  
  X_se = t(SpiecEasi::clr(counts_el_shuffled + 1, 1))
  graph.se = compute_graph_with_target_number_of_edges_mb(X_se, 1200)
  
  
  S_flashweave = flashweave_wrapper(counts_el_shuffled, 0.5)
  graph.flashweave = compute_graph_with_target_number_of_edges_corr(S_flashweave, target_edge_number)
  
  ########## Sparcc 
  
  scc = sparcc(counts_el_shuffled)
  graph.sparcc = compute_graph_with_target_number_of_edges_corr(abs(scc$Cor), target_edge_number)
  
  ######## MAGMA
  
  
  res_norm <- magma.norm(data=counts_el_shuffled, distrib="ZINB", seq_depth="GMPR")
  Z <- res_norm$z
  graph.magma = compute_graph_with_target_number_of_edges_mb(Z, 1200)
  
  
  z1 = graph.se != 0
  z2 = graph.magma != 0
  z3 = graph.zi != 0
  z4 = graph.sparcc != 0
  z5 = graph.flashweave != 0
  
  graph = list(z1, z2, z3, z4, z5)
  areas[k,] = draw_venn(graph, "venn_shuffled_counts.pdf")
}

##save(areas, file = "areas.rda")

intersections = areas
areas_mean = apply(areas, 2, mean)
i = 3
pdf(file="control_mean.pdf", height=11, width=11) 
plot.new()
venn.plot <- draw.quintuple.venn(
  area1 = intersections[i,1], 
  area2 = intersections[i,2], 
  area3 = intersections[i,3], 
  area4 = intersections[i,4], 
  area5 = intersections[i,5], 
  n12 = intersections[i,6], 
  n13 = intersections[i,7], 
  n14 = intersections[i,8], 
  n15 = intersections[i,9], 
  n23 = intersections[i,10],
  n24 = intersections[i,11],
  n25 = intersections[i,12],
  n34 = intersections[i,13],
  n35 = intersections[i,14],
  n45 = intersections[i,15],
  n123 = intersections[i,16],
  n124 = intersections[i,17],
  n125 = intersections[i,18],
  n134 = intersections[i,19],
  n135 = intersections[i,20],
  n145 = intersections[i,21], 
  n234 = intersections[i,22],
  n235 = intersections[i,23],
  n245 = intersections[i,24],
  n345 = intersections[i,25],
  n1234 = intersections[i,26],
  n1235 = intersections[i,27],
  n1245 = intersections[i,28],
  n1345 = intersections[i,29],
  n2345 = intersections[i,30],
  n12345 = intersections[i,31],
  category = c("Spiec-Easi", "Magma", "ZiLN", "Sparcc", "Flashweave"),
  fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"),
  cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), 
  cat.cex = 1.5,
  margin = 0.07,
  ind = TRUE
);
dev.off()


area1 = areas[i,1]
area2 = areas[i,2] 
area3 = areas[i,3] 
area4 = areas[i,4] 
area5 = areas[i,5] 
n12 = areas[i,6]
n13 = areas[i,7] 
n14 = areas[i,8] 
n15 = areas[i,9] 
n23 = areas[i,10]
n24 = areas[i,11]
n25 = areas[i,12]
n34 = areas[i,13]
n35 = areas[i,14]
n45 = areas[i,15]
n123 = areas[i,16]
n124 = areas[i,17]
n125 = areas[i,18]
n134 = areas[i,19]
n135 = areas[i,20]
n145 = areas[i,21] 
n234 = areas[i,22]
n235 = areas[i,23]
n245 = areas[i,24]
n345 = areas[i,25]
n1234 = areas[i,26]
n1235 = areas[i,27]
n1245 = areas[i,28]
n1345 = areas[i,29]
n2345 = areas[i,30]
n12345 = areas[i,31]


a31 <- n12345
a30 <- n1234 - a31
a29 <- n1235 - a31
a28 <- n1245 - a31
a27 <- n1345 - a31
a26 <- n2345 - a31
a25 <- n245 - a26 - a28 - a31
a24 <- n234 - a26 - a30 - a31
a23 <- n134 - a27 - a30 - a31
a22 <- n123 - a29 - a30 - a31
a21 <- n235 - a26 - a29 - a31
a20 <- n125 - a28 - a29 - a31
a19 <- n124 - a28 - a30 - a31
a18 <- n145 - a27 - a28 - a31
a17 <- n135 - a27 - a29 - a31
a16 <- n345 - a26 - a27 - a31
a15 <- n45 - a18 - a25 - a16 - a28 - a27 - a26 - a31
a14 <- n24 - a19 - a24 - a25 - a30 - a28 - a26 - a31
a13 <- n34 - a16 - a23 - a24 - a26 - a27 - a30 - a31
a12 <- n13 - a17 - a22 - a23 - a27 - a29 - a30 - a31
a11 <- n23 - a21 - a22 - a24 - a26 - a29 - a30 - a31
a10 <- n25 - a20 - a21 - a25 - a26 - a28 - a29 - a31
a9 <- n12 - a19 - a20 - a22 - a28 - a29 - a30 - a31
a8 <- n14 - a18 - a19 - a23 - a27 - a28 - a30 - a31
a7 <- n15 - a17 - a18 - a20 - a27 - a28 - a29 - a31
a6 <- n35 - a16 - a17 - a21 - a26 - a27 - a29 - a31
a5 <- area5 - a6 - a7 - a15 - a16 - a17 - a18 - a25 - 
  a26 - a27 - a28 - a31 - a20 - a29 - a21 - a10
a4 <- area4 - a13 - a14 - a15 - a16 - a23 - a24 - a25 - 
  a26 - a27 - a28 - a31 - a18 - a19 - a8 - a30
a3 <- area3 - a21 - a11 - a12 - a13 - a29 - a22 - a23 - 
  a24 - a30 - a31 - a26 - a27 - a16 - a6 - a17
a2 <- area2 - a9 - a10 - a19 - a20 - a21 - a11 - a28 - 
  a29 - a31 - a22 - a30 - a26 - a25 - a24 - a14
a1 <- area1 - a7 - a8 - a18 - a17 - a19 - a9 - a27 - 
  a28 - a31 - a20 - a30 - a29 - a22 - a23 - a12
areas <- c(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, 
           a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, 
           a21, a22, a23, a24, a25, a26, a27, a28, a29, a30, 
           a31)