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

###
# @author: Tobias Bauer
# Implementation of mass-action calculations based on
# section 6.2 "Effective abundance, pseudo-inverse and mass action" in
# the Supplementary Information of Staniczenko et al. 2013 Nature comm.
# "The ghost of nestedness in ecological networks"

massaction <- function(web) {
    # Computes the effective abundances for all species in web and
    # normalizes the matrix for mass action.
    # Args:
    #   web: the preference matrix as an integer matrix object.
    # Return:
    #  The preference matrix normalized for mass action.
    r <- nrow(web)
    c <- ncol(web)
    # Bipartite web as 1-D vector. Note: This is done by column i.e.
    # 30 0 16
    # 33 0  0 is converted to (30, 33, 19, 0, 0, 0, 1, 2, 16, 0, 15, 7)
    # 19 1 15
    #  0 2  7
    Bij <- as.vector(web)
    # logarithm of all non-zero vector elements
    y <- log(Bij[Bij != 0])
    # Create zero-filled matrix M
    M <- matrix(0, nrow = length(y), ncol = r + c)
    k <- 1
    # Fill matrix M
    # this has to be done in the same order as the web to 1-D vector conversion
    for (j in 1:c) {
        for (i in 1:r) {
            if (web[i, j] != 0) {
                M[k, i] <- 1
                M[k, r + j] <- 1
                # Increment k by one so the next row in M is filled
                k <- k + 1
            }
        }
    }
    # Solving the equation y = Mx for x
    # x = [x_i; x_j]
    # Use pseudo inverse (M^T*M)^(-1)*M^T * y = x
    x <- corpcor::pseudoinverse(M) %*% y
    # Applying exponential function to vector x,
    # since the logarithm was applied to y in the steps above
    x <- exp(x)
    # Devide each row in the incidence matrix by corresponding mass action in x.
    # MARGIN = 1 indicates applying FUN for each row in web.
    # The first r elements of x contain the mass action of "row"-species.
    gamma <- sweep(web, MARGIN = 1, STATS = x[1:r], FUN = "/")
    # Devide each column in the matrix by corresponding mass action in x.
    # The last c elements of x contain the mass action of "column"-species.
    gamma <- sweep(gamma, MARGIN = 2, STATS = x[(r + 1):(r + c)], FUN = "/")
    return(gamma)
}
