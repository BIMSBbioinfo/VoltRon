#' @include generics.R
#'
NULL

####
# Neighbor graphs ####
####

#' @param assay assay
#' @param method the method spatial connectivity
#' @param ... additional parameters passed to \code{FNN:get.knn}
#'
#' @rdname getSpatialNeighbors
#' @method getSpatialNeighbors VoltRon
#'
#' @importFrom tripack tri.mesh neighbours
#' @importFrom igraph add_edges simplify make_empty_graph vertices
#' @importFrom FNN get.knn
#'
#' @export
#'
getSpatialNeighbors.VoltRon <- function(object, assay = NULL, method = "delaunay", ...){

  # get coordinates
  coords <- vrCoordinates(object, assay = assay, reg = TRUE)

  # get spatial graph per assay
  assay_names <- vrAssayNames(object, assay = assay)
  spatialedges_list <- list()
  for(assy in assay_names){
    cur_coords <- coords[grepl(paste0(assy, "$"), rownames(coords)), ]
    spatialedges <-
      switch(method,
             delaunay = {
               tess <- tripack::tri.mesh(cur_coords[,1], cur_coords[,2])
               nnedges <- tripack::neighbours(tess)
               nnedges <- mapply(function(x,y) {
                 do.call(c,lapply(y, function(z) c(x,z)))
               }, 1:length(nnedges), nnedges)
               nnedges <- unlist(nnedges)
               nnedges <- rownames(cur_coords)[nnedges]
               nnedges
             },
             NN = {
               nnedges <- FNN::get.knn(cur_coords, k = 1)
               nnedges <- cbind(1:nrow(cur_coords), nnedges)
               nnedges <- apply(nnedges, 1, function(x){
                 do.call(c,lapply(x[-1], function(y) return(c(x[1],y))))
               })
               nnedges <- unlist(nnedges)
               nnedges <- rownames(cur_coords)[nnedges]
               nnedges
             })
    spatialedges_list[[assy]] <- spatialedges
  }
  spatialedges <- unlist(spatialedges_list)

  # make graph and add edges
  graph <- make_empty_graph(directed = FALSE) + vertices(rownames(coords))
  graph <- add_edges(graph, edges = spatialedges)
  graph <- simplify(graph, remove.multiple = TRUE, remove.loops = FALSE)
  vrGraph(object, assay = assay, graph.type = "delaunay") <- graph

  # return
  return(object)
}

####
# Neighbor test ####
####

vrNeighbourhoodEnrichment <- function(object, assay = NULL, assay.type = NULL, group.by = NULL, method = "dixon", graph.type = "delaunay", seed = 1, num.sim = 100){

  # set the seed
  set.seed(seed)

  # check object
  if(!inherits(object, "VoltRon"))
    stop("Please provide a VoltRon object!")

  # sample metadata
  sample.metadata <- SampleMetadata(object)

  # get assay names
  assay_names <- vrAssayNames(object, assay = assay)

  # get entity type and metadata
  if(is.null(assay.type)){
    assay_types <- vrAssayTypes(object, assay = assay)
    if(length(unique(assay_types)) == 1){
      metadata <- Metadata(object, type = unique(assay_types))
    } else {
      stop("Please select assay.type as 'cell', 'spot' or 'ROI'")
    }
  } else {
    metadata <- Metadata(object, type = assay.type)
  }

  # coords
  coords <- as.data.frame(vrCoordinates(object))

  dixon_results <- NULL
  if(method == "dixon"){
    datax <- cbind(coords, cur_metadata[[group.by]])
    dixon_results <- dixon::dixon(datax)
    return(dixon_results$tablaZ)
  } else if(method == "permutation"){
    graph <- vrGraph(object, graph.type = graph.type)
    neighbors_graph <- igraph::neighborhood(graph)


  }

  return(dixon_results)
}

