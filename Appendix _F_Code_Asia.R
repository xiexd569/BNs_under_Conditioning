# Load necessary libraries
library(bnlearn)
library(igraph)

# Function to get ancestors including the node itself
ancestors_with_self <- function(G, node) {
  # Get ancestors of the node
  ancestors_set <- bnlearn::ancestors(G, node)
  # Include the node itself in the ancestor set
  return(union(ancestors_set, node))
}

# Function to find minimal fill edge set
minimal_fill_edge_set <- function(G, alpha, C) {
  # Initialize
  topological_order <- alpha
  nodes <- C 
  E_0 <- vector("list", length = 0)  # Initialize empty list
  
  # Update G by removing children of nodes in C
  for (v in C) {
    ch <- bnlearn::children(G, v)
    for (u in ch) {
      G <- drop.arc(G, v, u)
    }
  }
  
  # Get ancestors of nodes in C
  ancestors_list <- lapply(nodes, function(node) ancestors_with_self(G, node))
  ancestor_union_all <- Reduce(union, ancestors_list)
  nodes_of_interest <- ancestor_union_all 
  G_0 <- induced_subgraph(as.igraph(G), nodes_of_interest)  # Induced subgraph
  G_00 <- as.bn(G_0)
  
  # Order vertices according to alpha
  k <- length(nodes(G_00))
  vertices <- intersect(rev(topological_order), nodes_of_interest)
  v_structures <- list()
  
  # If number of vertices >= 3, perform fill algorithm
  if (k >= 3) {
    for (v in vertices[1:(k-2)]) {
      parents <- bnlearn::parents(G_00, v)
      if (length(parents) >= 2) {
        comb <- combn(parents, 2)
        for (i in 1:ncol(comb)) {
          x <- comb[1, i]
          y <- comb[2, i]
          if (all(!(x %in% bnlearn::parents(G_00, y)) & !(y %in% bnlearn::parents(G_00, x)))) {
            if (which(topological_order == x) < which(topological_order == y)) {
              G_00 <- set.arc(G_00, x, y, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
              edge1 <- paste(x, y, sep = "->")
              E_0[[length(E_0) + 1]] <- edge1    # Add element
            } else {
              G_00 <- set.arc(G_00, y, x, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
              edge1 <- paste(y, x, sep = "->")
              E_0[[length(E_0) + 1]] <- edge1    # Add element
            }
          }
        }
      }
    }
  }  
  if (length(E_0) != 0) {
    E_0 <- t(unlist(E_0))                 # Convert to vector 
  }
  
  for (v in C) {
    G_00 <- remove.node(G_00, v)
  }
  
  # Return results
  list(E_alpha_R = E_0, G_alpha_R = G_00)
}

# Function to find minimal fill edge set using l_m_p algorithm
minimal_fill_edge_set_l_m_p <- function(G, alpha, C) {
  dag <- G
  topological_order <- alpha
  nodes <- C
  E_0 <- vector("list", length = 0)  # Initialize empty list
  
  # Get ancestors of nodes in C
  ancestors_list <- lapply(nodes, function(node) ancestors_with_self(G, node))
  ancestor_union_all <- Reduce(union, ancestors_list)
  nodes_of_interest <- ancestor_union_all 
  G_0 <- induced_subgraph(as.igraph(G), nodes_of_interest)  # Induced subgraph
  G_00 <- as.bn(G_0)
  
  # Order vertices according to alpha
  vertices <- intersect(rev(topological_order), nodes_of_interest)
  vertices3 <- vertices
  
  if (length(nodes_of_interest) >= 2) {
    for (v in vertices) {
      vertices1 <- vertices3
      parents1 <- bnlearn::parents(G_00, v)
      children1 <- bnlearn::children(G_00, v)
      p_c1 <- union(parents1, children1)
      mb1 <- bnlearn::mb(G_00, v)
      p_c_c1 <- union(p_c1, C)
      D <- setdiff(mb1, p_c_c1)
      
      if (length(D) > 0) {
        all_combinations <- c(list(NULL), unlist(lapply(1:length(D), function(x) combn(D, x, simplify = FALSE)), recursive = FALSE))
      } else {
        vertices3 <- setdiff(vertices3, v)
        next
      }
      
      for (c_set in all_combinations) {
        if (length(c_set) > 0) {
          s <- union(c_set, parents1)
          c1 <- union(s, C)
          vertices2 <- setdiff(vertices1, c(v, s))
          
          if (length(vertices2) > 0) {
            d_tf3 <- sapply(vertices2, function(z) bnlearn::dsep(G_00, v, z, c1))
          } else {
            d_tf3 <- TRUE
          }
          
          if (all(d_tf3)) {
            for (z in c_set) {
              G_00 <- set.arc(G_00, z, v, check.cycles = TRUE, check.illegal = TRUE)
              edge1 <- paste(z, v, sep = "->")
              E_0[[length(E_0) + 1]] <- edge1    # Add element
            }
          }
        }
      } 
      vertices3 <- setdiff(vertices3, v) 
    }
  }
  
  if (length(E_0) != 0) {
    E_0 <- t(unlist(E_0))                 # Convert to vector 
  }
  
  for (v in C) {
    G_00 <- remove.node(G_00, v)
  }
  
  list(E_alpha_R = E_0, G_alpha_R = G_00)
}

# Load data
load("D:/Users/xie_x/Downloads/Asia.rda")

# Define subset C
C <- c("xray", "dysp")

# Convert network to igraph
network_string <- modelstring(bn)
asianet <- model2network(network_string)
asianet1 <- as.igraph(asianet)

# Set seed for reproducibility
set.seed(1232024)

# Get topological order
alpha <- node.ordering(asianet, debug = FALSE)

# Measure execution time for minimal_fill_edge_set
start_time1 <- Sys.time()
for (i in 1:100) {
  z1 <- minimal_fill_edge_set(asianet, alpha, C)
}
end_time1 <- Sys.time()
execution_time1 <- as.numeric(difftime(end_time1, start_time1, units = "secs")) / 100

# Measure execution time for minimal_fill_edge_set_l_m_p
start_time1_d_1 <- Sys.time()
for (i in 1:100) {
  z1_d_1 <- minimal_fill_edge_set_l_m_p(asianet, alpha, C)
}
end_time1_d_1 <- Sys.time()
execution_time1_d_1 <- as.numeric(difftime(end_time1_d_1, start_time1_d_1, units = "secs")) / 100

# Print execution times
print(execution_time1)
print(execution_time1_d_1)

# Compare results
print(setdiff(z1$E_alpha_R, z1_d_1$E_alpha_R))
print(setdiff(z1_d_1$E_alpha_R, z1$E_alpha_R))

# Plot the graph with different edges highlighted
G <- asianet
nodes <- C
ancestors_list <- lapply(nodes, function(node) ancestors_with_self(G, node))
ancestor_union_all <- Reduce(union, ancestors_list)
nodes_of_interest <- ancestor_union_all 

G_0 <- induced_subgraph(as.igraph(G), nodes_of_interest) 
g1 <- G_0
g2 <- as.igraph(z1$G_alpha_R)

edges_g1 <- E(g1)
edges_g2 <- E(g2)

diff_edges <- list()
for (edge in edges_g2) {
  if (!(edge %in% edges_g1)) {
    diff_edges <- append(diff_edges, edge)
  }
}

diff_edges_indices <- which(edges_g2 %in% diff_edges)

edge_colors <- rep("black", ecount(g2))  # Default color
edge_linetypes <- rep("solid", ecount(g2))  # Default linetype

edge_colors[diff_edges_indices] <- "red"
edge_linetypes[diff_edges_indices] <- "dashed"

layout_matrix <- layout_with_fr(g2, niter = 10000)

plot(g2, layout = layout_matrix, 
     vertex.size = 30,  # Increase node size
     vertex.label.cex = 1.0,  # Adjust label size
     edge.arrow.size = 0.3,  # Adjust arrow size
     edge.color = edge_colors, edge.lty = edge_linetypes, main = "DAG 2 with Different Edges Highlighted")