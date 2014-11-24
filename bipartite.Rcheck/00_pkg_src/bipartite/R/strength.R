strength <- function(web, type="Bascompte"){
	# there are two definitions of STRENGTH, that of Bascompte et al. 2006 as the sum of dependencies of a species, and that of Barrat et al. 2004 as the weighted sum of links.
    # Bascompte et al.'s strength sums to the number of species in the other group
    # Barrat's strength is simply the number of interactions, a trivial measure of a species importance; in contrast to the claim of Poisot et al. (2012, p. 1155), this definition of strength gives no information of the centrality of a species within a network structure (and neither does the other).
    # Default: Bascompte
	if (!(type %in% c("Bascompte", "Barrat"))) stop("Argument 'type' not recognised: should be either 'Bascompte' or 'Barrat'.")

	if (type == "Bascompte"){
		depL <- web/matrix(rowSums(web), nrow=NROW(web), ncol=NCOL(web), byrow=FALSE)
		SH <- colSums(depL)
	}
		
	if (type == "Barrat"){
		SH <- colSums(web)/sum(web)
	}
	
	return(SH)
}
#data(Safariland)
#s1 <- strength(Safariland, type="Barrat")
#s2 <- strength(Safariland, type="Bascompte")
#plot(s1, s2, log="x")
#cor.test(s1, s2, type="ken")