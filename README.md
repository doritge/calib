# Cellular Automata
Generating artwork using various Cellular Automata and Elementary Cellular Automata rules

## General
Cellular Automata (CA) are mathematical models of computation that transform a grid of cells, each in a specific state, over discrete time steps. The transformation is based on rules that depend only on the states of the cell and its neighbors.

The code is implemented in R, sometimes using C++ code for enhanced performance

## Two-dimensional rules
Scripts cell_automata.R and ca_iterators.cpp

**Totalistic CA**

A Totalistic CA rule defines the new state of a cell based on the sum/average/count of the state of its neighboring cells. This implementation uses the sum of the states of the cell and its neighboring cells as an index into a look up table of possible states.

<img src="/examples/total6242M.png" width="400" height="400" />

**Index CA**

An Index CA rule uses the state of the current cell as an index to a list of possible new states. In this implementation the state of the cell serves as an index to the list of states of its current neighbors.

<img src="/examples/index7633.png" width="400" height="400" />

**Cyclic CA**

The Cyclic CA rule increments the state of the current cell when the count of the neighboring cells that are exactly 1 above the state of the current cell exceeds a threshold.

<img src="/examples/cyclic4453.png" width="400" height="400" />

## One-dimensional rules
Scripts elementary.R

**Wolfram binary**

The original Wolframs's 1-dimensional, 2-state CA. The 3 binary states of the cell and its two neighbors form a vector of 8 possible binary values, which in turn result in 8 binary results, i.e. 256 possible rules. The rule number is defined by the binary representation of the states outcome in each such vector. Once a rule has been selected, the index into this outcome vector is computed as binary representation of the current state of the cell and its two surrounding neighbors.

<img src="/examples/wolfram30.png" width="400" height="400" />

**Wolfram totalistic**

The sum of the cell and its 2 neighbors serves as the index to a vector of new states.

<img src="/examples/Etotal7648R.png" width="400" height="400" />

**Elementary muli-state**

The multi-state cell is mapped to a binary cell and the next binary row is calculated by using the binary rule as in the *Wolfram binary* rule. The new binary row is now used as control for generating the new multi-state row from the current multi-state row. If the new cell is 0, it is copied as is. Otherwise, in the "S" method, the current cell is copied to the new cell and in the "N" method, the next state is sampled out of the current cell and its neighbors.

<img src="/examples/multi73_9619.png" width="400" height="400" />

