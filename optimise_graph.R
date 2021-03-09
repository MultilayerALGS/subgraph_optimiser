library(reticulate)

# W_filename and data_filename are actual names to files.
optimise_graph_from_file <- function(W_filename, data_filename) {
  utils <- import_from_path("utils")
  optimiser <- import_from_path("optimiser")
  data <- utils$read_data(data_filename)
  graph <- utils$read_matrix(W_filename)
  opt_object <- optimiser$Optimiser(graph, data)
  optimal_graph <- opt_object$iterative_opt(remove=TRUE, add=FALSE)
  optimal_graph
}
# The result from running this function is a "Graph" object from Python. It has
# some usual graph functions, the one that's probably relevant is graph$edges()
# which gives a list of the edges still ini the graph.
# > optimal_graph$edges()


# This function uses data from R, rather than reading from files.
# W_matrix needs to be either a list of lists of zeros and ones,
# or it needs to be a matrix of zeros and ones.
# data_list needs to be a list of numbers.
optimise_graph_from_data <- function(W_matrix, data_list) {
  if (all(class(W_matrix) == c("matrix", "array"))) {
    W_matrix = lapply(apply(W_matrix, 1, as.list), unlist) 
  }
  optimiser <- import_from_path("optimiser")
  opt_object <- optimiser$Optimiser(W_matrix, data_list)
  optimal_graph <- opt_object$iterative_opt(remove=TRUE, add=FALSE)
  optimal_graph
}

W <- list(list(0,1,0,0,0,0),
          list(1,0,1,0,0,0),
          list(0,1,0,1,0,0),
          list(0,0,1,0,1,0),
          list(0,0,0,1,0,1),
          list(0,0,0,0,1,0))

data <- list(100, 99, 0, 100, 199, 200)
optimise_graph_from_data(W, data)
