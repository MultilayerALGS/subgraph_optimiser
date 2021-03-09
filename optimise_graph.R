library(reticulate)

# W_filename and data_filename are actual names to files.
optimise_graph <- function(W_filename, data_filename) {
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