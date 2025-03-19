##Finding ancestors and themselves
ancestors_with_self <- function(G, node) {
  # Retrieve the ancestors of the given node
  ancestors_set <- bnlearn::ancestors(G, node)
  
  # Include the node itself in the ancestor set
  return(union(ancestors_set, node))
}

##Finding minimal filling edge set by the ROS algorithm
minimal_fill_edge_set <- function(G, alpha, C) {
  
  topological_order <- alpha
  nodes <- C 
  E_0 <- vector("list", length = 0)  # Initialize an empty list
  
  # Update graph by removing outgoing edges from nodes in C
  for (v in C) {
    ch <- bnlearn::children(G, v)
    for (u in ch) {
      G <- drop.arc(G, v, u)
    }
  }
  
  # Compute the set of ancestors for the nodes in C
  ancestors_list <- lapply(nodes, function(node) ancestors_with_self(G, node))
  ancestor_union_all <- Reduce(union, ancestors_list)
  nodes_of_interest <- ancestor_union_all 
  
  # Induce a subgraph containing only the relevant nodes
  G_0 <- induced_subgraph(as.igraph(G), nodes_of_interest)  
  G_00 <- as.bn(G_0)
  
  # Order vertices according to the given topological order
  k <- length(nodes(G_00))
  vertices <- intersect(rev(topological_order), nodes_of_interest)
  v_structures <- list()
  
  # Apply minimal fill-in algorithm if there are at least 3 vertices
  if (k >= 3) {
    for (v in vertices[1:(k-2)]) {
      parents <- bnlearn::parents(G_00, v)
      if (length(parents) >= 2) {
        comb <- combn(parents, 2)
        for (i in 1:ncol(comb)) {
          x <- comb[1, i]
          y <- comb[2, i]
          
          # Ensure that x and y are not already connected
          if (all(!(x %in% bnlearn::parents(G_00, y)) & !(y %in% bnlearn::parents(G_00, x)))) {
            if (which(topological_order == x) < which(topological_order == y)) {
              G_00 <- set.arc(G_00, x, y, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
              edge1 <- paste(x, y, sep = "->")
              E_0[[length(E_0) + 1]] <- edge1  # Add edge to list
            } else {
              G_00 <- set.arc(G_00, y, x, check.cycles = TRUE, check.illegal = TRUE, debug = FALSE)
              edge1 <- paste(y, x, sep = "->")
              E_0[[length(E_0) + 1]] <- edge1  # Add edge to list
            }
          }
        }
      }
    }
  }  
  
  # Convert edge list to vector if not empty
  if (length(E_0) != 0) {
    E_0 <- t(unlist(E_0))
  }
  
  # Remove nodes in C from the final graph
  for (v in C) {
    G_00 <- remove.node(G_00, v)
  }
  
  # Return the modified graph and the added edges
  list(E_alpha_R = E_0, G_alpha_R = G_00)
}
