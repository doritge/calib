##====================================
## Elementary Cellular Automata
##====================================

# calculate left & right neighborhoods with wrapping
left_neighbors <- function(vect){
  n <- length(vect)
  left <- dplyr::lag(vect)
  left[1] <- vect[n]
  return(left)
}
right_neighbors <- function(vect){
  n <- length(vect)
  right <- dplyr::lead(vect)
  right[n] <- vect[1]
  return(right)
}

# convert rule number to binary rule array
rule2vect <- function(rule_num){
  rule <- c()
  for (i in 1:8){
    rule <- c(rule, rule_num%%2)
    rule_num <- rule_num%/%2
  }
  return(rule)
}

#' @title Wolfram CA binary
#' @description Generate the basic binary Wolfram Cellular Automata
#'
#' @param rule_num Wolfram rule number
#' @param ncolumns number of columns in grid
#' @param nrows number of rows in grid
#' @param first_mode mode of generation of first row. ("M": mid row on, "R": random)
#' @param seed seed, defaults to NULL
#'
#' @return A matrix of cells
#'
#' @export
#'

wolfram <- function(rule_num,
                    ncolumns, nrows,
                    first_mode,
                    seed = NULL){

  # calculate next row
  next_row <- function(row, rule){
    row_left <- left_neighbors(row)
    row_right <- right_neighbors(row)
    indices <- row_left + 2*row + 4*row_right + 1
    new_wolfram <- rule[indices]
  }

  if(!is.null(seed)) set.seed(seed)
  rule_vect <- rule2vect(rule_num)

  # generate first row according to mode
  if (first_mode == "M"){
    first_row <- rep(0, ncolumns)
    first_row[ncolumns%/%2] <- 1
  }
  if (first_mode == "R"){
    first_row = sample(0:1, ncolumns, replace = TRUE)
  }

  # iterate on rows to create a 2-dimensional cell grid
  mat <- purrr::accumulate(1:(nrows-1), ~next_row(.x, rule = rule_vect), .init = first_row) |>
    purrr::reduce(rbind)

  return(mat)
}

#' @title Wolfram CA totalistic
#' @description Generate the totalistic Wolfram Cellular Automata. It uses summation rule and it supports multiple states
#'
#' @param rule_num if NULL - random generation of rule, otherwise generated based on # of states
#' @param ncolumns number of columns in grid
#' @param nrows number of rows in grid
#' @param nstates number of states
#' @param first_mode mode of generation of first row. ("M": mid row on, "R": random)
#' @param seed seed, defaults to NULL
#'
#' @return A matrix of cells
#'
#' @export
#'

wolfram_total <- function(rule_num = NULL,
                          ncolumns, nrows,
                          nstates,
                          first_mode,
                          seed = NULL){

  next_row <- function(row, rule, ns){
    sum <- left_neighbors(row) + row + right_neighbors(row)
    new_total <- rule[sum + 1]
    return(new_total)
  }

  # initialize parameters
  if(!is.null(seed)) set.seed(seed)

  ## generate rule vector
  if(is.null(rule_num)){
    rule_vect <- sample(0:(nstates-1), 3*nstates - 2, replace = TRUE)     # random
  }else{
    rule_vect <- c()
    for (i in 1:(3*nstates - 2)){
      rule_vect <- c(rule_vect, rule_num%%nstates)
      rule_num <- rule_num%/%nstates
    }
  }

  ## generate first row based on mode. M: mid cell, R: random
  if (first_mode == "M"){
    first_row <- rep(0, ncolumns)
    first_row[ncolumns%/%2] <- 1
  }
  if (first_mode == "R"){
    first_row = sample(0:(nstates-1), ncolumns, replace = TRUE)
  }

  # iterate on rows to create a 2-dimensional cell grid
  mat <- purrr::accumulate(1:(nrows-1), ~next_row(.x, rule = rule_vect),
                           .init = first_row) |>
    purrr::reduce(rbind)

  return(mat)
}

#' @title Elementary CA Multi-state
#' @description Generate a multi-state Elementary Cellular Automata by mapping the multi-state cells to binary cells that are used as controls for the state selection
#'
#' @param rule_num Wolfram rule number
#' @param ncolumns number of columns in grid
#' @param nrows number of rows in grid
#' @param nstates number of states
#' @param first_mode mode of generation of first row. ("M": mid row on, "R": random)
#' @param gen_mode mode of state selection. ("S": state of current cell, "N": sample of cell state and neighbors')
#' @param seed seed, defaults to NULL
#'
#' @return A matrix of cells
#'
#' @export
#'

elementary_multi <- function(rule_num,
                             ncolumns, nrows,
                             nstates,
                             first_mode,
                             gen_mode,
                             seed = NULL){

  # calculate next row using Wolfram 2-state rule as control
  next_row <- function(row, rule){
    # calculate row neighborhoods
    row_left <- left_neighbors(row)
    row_right <- right_neighbors(row)

    # calculate Wolfram and its neighborhoods
    binary <- ifelse(row == 0, 0, 1)      # map row to binary
    binary_left <- left_neighbors(binary)
    binary_right <- right_neighbors(binary)

    # generate a new binary row
    indices <- binary_right + 2*binary + 4*binary_left + 1
    new_binary <- rule[indices]

    # generate a new row
    new_row <- rep(0, ncolumns)
    for (i in 1: ncolumns){
      if (new_binary[i] != 0){       # new binary cell is alive, now generate its state

        if (gen_mode == "S"){
          # consider only current cell
          if (row[i] == 0) new_row[i] <- sample(1:(nstates-1), 1)
          else new_row[i] <- row[i]
        }

        if (gen_mode == "N"){
          # consider current cell & 2 neighbors
          states_buffer <- c()
          if (binary_left[i] != 0) states_buffer <- c(states_buffer, row_left[i])
          if (binary[i] != 0) states_buffer <- c(states_buffer, row[i])
          if(binary_right[i] !=0) states_buffer <- c(states_buffer, row_right[i])
          if (is.null(states_buffer)) new_row[i] <- sample(0:(nstates-1), 1)
          else new_row[i] <- sample(states_buffer, 1)
        }

      }
    }
    return(new_row)
  }

  if(!is.null(seed)) set.seed(seed)
  rule_vect <- rule2vect(rule_num)

  ## generate first row based on mode. M: mid cell, R: random
  if (first_mode == "M"){
    first_row <- rep(0, ncolumns)
    first_row[ncolumns%/%2] <- 1
  }
  if (first_mode == "R"){
    first_row = sample(0:(nstates-1), ncolumns, replace = TRUE)
  }

  # iterate on rows to create a 2-dimensional cell matrix
  mat <- purrr::accumulate(1:(nrows-1), ~next_row(.x, rule = rule_vect),
                           .init = first_row) |>
    purrr::reduce(rbind)

  return(mat)
}
