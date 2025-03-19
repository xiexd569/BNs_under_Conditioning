##Finding ancestors and themselves
ancestors_with_self <- function(G, node) {
  # Retrieve the ancestors of the given node
  ancestors_set <- bnlearn::ancestors(G, node)
  
  # Include the node itself in the ancestor set
  return(union(ancestors_set, node))
}


##finding minimal filling edge set by the local Markov property
minimal_fill_edge_set_l_m_p <- function(G, alpha, C) {
  
  dag <- G
  topological_order <- alpha
  nodes <- C
  E_0 <- vector("list", length = 0)  # Initialize an empty list
  
  # Compute the set of ancestors for nodes in C
  ancestors_list <- lapply(nodes, function(node) ancestors_with_self(G, node))
  ancestor_union_all <- Reduce(union, ancestors_list)
  nodes_of_interest <- ancestor_union_all 
  
  # Induce a subgraph containing only the relevant nodes
  G_0 <- induced_subgraph(as.igraph(G), nodes_of_interest)
  G_00 <- as.bn(G_0)
  
  # Order vertices according to the given topological order
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
      
      # Iterate over all subsets of D
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
              E_0[[length(E_0) + 1]] <- edge1  # Add edge to list
            }
          }
        }
      } 
      vertices3 <- setdiff(vertices3, v) 
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