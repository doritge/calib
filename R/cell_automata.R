##====================================
## Cellular Automata with various rules
##====================================

#' @title Totalistic Cellular Automata
#' @description Generate a Totalistic CA with animation option
#'
#' @param ncolumns number of columns in grid
#' @param nrows number of rows in grid
#' @param init_mode mode of initial configuration, "R": random, "M": middle of grid
#' @param neighborhood configuration of neighborhood (x, y, weight)
#' @param iterations number of iterations
#' @param animate animate iterations? defaults to FALSE
#' @param seed seed, defaults to NULL
#'
#' @return A a single long grid matrix or frames of long grid for animation
#'
#' @export
#'

ca_total <- function(ncolumns, nrows,
                     nstates,
                     init_mode,
                     neighborhood,
                     iterations,
                     animate = FALSE,
                     seed = NULL){

  if(!is.null(seed)) set.seed(seed)

  # create next generation look up table
  rule_vect <- sample(0:(nstates-1), nrow(neighborhood)*nstates-1, replace = TRUE)

  ## generate initial matrix
  if(init_mode == "R"){
    mat <- ambient::long_grid(x = seq(0, 1, length = ncolumns), y = seq(0, 1, length = nrows)) |>
      dplyr::mutate(v = sample(0:(nstates-1), ncolumns*nrows, replace = TRUE)) |>
      long2mat()
  }
  if(init_mode == "M"){
    mat <- matrix(rep(0, ncolumns*nrows), nrow = nrows)
    mat[nrows/2, ncolumns/2] <- 1
  }

  # iterate on matrix
  if (animate == TRUE){
    anim <- mat2long(mat) |>
      mutate(frame = 0)
    for (f in 1:iterations){
      mat <- iterate_total(mat, neighborhood, rule_vect)
      long <- mat2long(mat) |>
        mutate(frame = f)
      anim <- rbind(anim, long)
    }
    return(anim)
  }else{
    for(i in 1:iterations){
      mat <- iterate_total(mat, neighborhood, rule_vect)
    }
    return(mat2long(mat))
  }
}

#' @title Cyclic Cellular Automata
#' @description if number of neighbors with state one unit greater than current is larger than the threshold increment current state by one
#' @param ncolumns number of columns in grid
#' @param nrows number of rows in grid
#' @param nstates number of states
#' @param neighborhood configuration of neighborhood (x, y, weight)
#' @param th threshold, should be proportionate to probability of each state
#' @param iterations number of iterations
#' @param animate animate iterations? defaults to FALSE
#' @param seed seed, defaults to NULL
#'
#' @return A a single long grid matrix or frames of long grid for animation
#'
#' @export
#'

ca_cyclic <- function(ncolumns, nrows,
                      nstates,
                      neighborhood,
                      th,
                      iterations,
                      animate = FALSE,
                      seed = NULL){

  if(!is.null(seed)) set.seed(seed)

  # create next generation lookup table
  rule_vect <- sample(0:(nstates-1), nrow(neighborhood)*nstates-1, replace = TRUE)

  ## generate initial matrix
  mat <- ambient::long_grid(x = seq(0, 1, length = ncolumns), y = seq(0, 1, length = nrows)) |>
    dplyr::mutate(v = sample(0:(nstates-1), ncolumns*nrows, replace = TRUE)) |>
    long2mat()

  # iterate on matrix
  if (animate == TRUE){
    anim <- mat2long(mat) |>
      mutate(frame = 0)
    for (f in 1:iterations){
      mat <- iterate_cyclic(mat, neighborhood, states = nstates, threshold = th)
      long <- mat2long(mat) |>
        mutate(frame = f)
      anim <- rbind(anim, long)
    }
    return(anim)
  }else{
    for(i in 1:iterations){
      mat <- iterate_cyclic(mat, neighborhood, states = nstates, threshold = th)
    }
    return(mat2long(mat))
  }
}

#' @title Index Cellular Automata
#' @description the state of the cell acts as an index to neighbor from which the cell's next state is taken.
#' @param nrows number of rows in grid
#' @param nstates number of states
#' @param neighborhood configuration of neighborhood (x, y, weight)
#' @param iterations number of iterations
#' @param seed seed, defaults to NULL
#'
#' @return A a single long grid matrix or frames of long grid for animation
#'
#' @export
#'

ca_index <- function(ncolumns, nrows,
                     nstates,
                     neighborhood,
                     iterations,
                     seed = NULL){

  if(!is.null(seed)) set.seed(seed)

  # create next generation lookup table
  rule_vect <- sample(0:(nrow(neighborhood)-1), nstates)

  ## generate initial matrix
  mat <- ambient::long_grid(x = seq(0, 1, length = ncolumns), y = seq(0, 1, length = nrows)) |>
    dplyr::mutate(v = sample(0:(nstates-1), ncolumns*nrows, replace = TRUE)) |>
    long2mat()

  # iterate on matrix
  for(i in 1:iterations){
    mat <- iterate_index(mat, neighborhood, R = rule_vect)
  }

  return(mat2long(mat))
}
