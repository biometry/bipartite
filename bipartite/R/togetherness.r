togetherness <- function(web, normalise=TRUE, FUN=mean, ...){
    # calculates the T-score="togetherness" for all pollinator species; the T-score represents
    # the average number of specis pair identical co-occurrences and co-absences.
    # (Stone & Roberts 1992)
    # J*(N-J-(A-J)-(B-J)) = J*(N-A-B+J)
    # for each species pair, we count the number of island pairs of the pattern (0,0,1,1) or (1,1,0,0,)
    # ... to be passed on to FUN (If a species occurs on all sites, then its distance
    # Carsten F. Dormann, Jan. 2008
    # didn't find a simple way to calculate the upper limit for togetherness.
    web <- web>0 # this whole concept works only on binary data!
    D <- designdist(t(web), method="J*(P-A-B+J)", terms="minimum", name="togetherness")
    # The minimum value for Ds is 0, for the special case were all species use the
    # hosts exactly co-occurringly.
    # The maximum value possible for each species is simply the product of number of 
    # 0s and 1s, so the maximum of a species pair is the minimum of each species maximum.
    if (normalise){
      maxD <- designdist(t(web), method="min(max(A*(P-A)), max(B*(P-B)))", terms="minimum")
      D <- D/maxD
    }
    
    FUN(D, ...)
}
# example:                                                                            
#m <- matrix(c(1,0,0, 1,1,0, 1,1,0, 0,1,1, 0,0,1), 5,3,TRUE)
#togetherness(m, normalise=FALSE)
#togetherness(m, normalise=TRUE)