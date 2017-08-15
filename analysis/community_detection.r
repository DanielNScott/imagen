detect_communities <- function(dist) {

  library(igraph)

  # Create graph from thresholded distances, remove multi-edges and loops
  g <- graph_from_adjacency_matrix(dist < 25, 'upper')
  g <- simplify(g)

  # Newman-Girvan algorithm
  ebc <- edge.betweenness.community(g, directed=F)

  # Now we have the merges/splits and we need to calculate the modularity
  # for each mergei.For this we'll use a function that for each edge
  # removed will create a second graph, check for its membership and use
  # that membership to calculate the modularity
  #mods <- sapply(0:ecount(g), function(i) {
  #  g2 <- delete.edges(g, ebc$removed.edges[seq(length=i)])
  #  cl <- clusters(g2)$membership

    # Compute modularity on the original graph g
  #  return(modularity(g,cl))
  #})

  # Now, let's color the nodes according to their membership
  #g2 <- delete.edges(g, ebc$removed.edges[seq(length=which.max(mods)-1)])
  #V(g)$color=clusters(g2)$membership

  # Let's choose a layout for the graph
  #g$layout <- layout.fruchterman.reingold

  #return()
}

# we can now plot all modularities
#plot(mods, pch=20)

# plot it
#plot(g, vertex.label=NA)

cov_dist <- function(A, B) {
  # Returns the Forstner-Moonen 1999 distance between covariance
  # matrices.
  #
  # Args:
  #   A, B covariance matrices
  #
  # Requires library 'geigen'

  library('geigen')

  eigen_vals <- geigen(A, B, symmetric = TRUE, only.values = TRUE)$values
  distance   <- sqrt(sum( sapply(eigen_vals, function(x){log(x)^2}) ))

  return(distance)
}
