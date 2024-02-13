###
# @author: Tobias Bauer
# Implementation of spectral radius as a nestedness measurement
# described in Staniczenko et al. 2013 Nature comm.
# "The ghost of nestedness in ecological networks"

spectral.radius <- function(web, mass.action.norm = FALSE){
    # If specified normalize the matrix for mass action
    if (mass.action.norm){
        web <- massaction(web)
    }
    # Convert the input bipartite web to a one-mode matrix
    net <- as.one.mode(web)
    # Calculate the Eigenvalues of the one-mode adjacency matrix
    ev <- eigen(net)$values
    # The spectral radius is the largest Eigenvalue,
    # which is the first vector element as they are
    # in descending order.
    return(ev[1])
}