webs2array <- function(...){
	# function to put several webs into an array, blowing the dimensions up to the union of species
	# makes sense only for webs featuring overlapping communities!
	webinput <- substitute(list(...))
	web.names <- as.character(webinput)[-1L]
	webs <- eval(webinput)
  if (is(webs[[1]], "list")){
    webs <- webs[[1]]
    web.names <- names(webs)
    if (is.null(web.names)) {web.names <- paste0("web",1:length(webs))}
  }
	
	if (length(webs) < 2) stop("Please provide more than one web as input!")
	
	all.names.HL <- sort(unique(unlist(lapply(webs, colnames)))) # this was sapply before, which caused duplicate names (unique didn't work)
	all.names.LL <- sort(unique(unlist(lapply(webs, rownames)))) # was sapply before
	
	web.array <- array(0, dim=c(length(all.names.LL), length(all.names.HL), length(web.names)), dimnames=list(all.names.LL, all.names.HL, web.names))
	for (i in 1:length(web.names)){
		ri <- rownames(webs[[i]]) # which was used before here, causing wrong ordering!
		ci <- colnames(webs[[i]])
		web.array[ri, ci, i] <- webs[[i]]
	}
	
	return(web.array)
}

#data(Safariland, vazquenc, vazquec)
#allin1 <- webs2array()
#allin1 <- webs2array(Safariland)
#allin1 <- webs2array(Safariland, vazquenc, vazquec)
## now we can compute distance between two webs:
#vegdist(t(cbind(as.vector(allin1[,,2]), as.vector(allin1[,,3]))), method="jacc")

# new option: input is a weblist
# allin1 <- webs2array(list(Safariland, vazquenc)) # list has no names: new names "web1" etc. are given for third dimension of webarray
# allin1 <- webs2array(list(Safariland=Safariland, Vazquenc=vazquenc)) # if present, list names are used

## to get the input part sorted out:
#webs2array <- function(...){
#	out <- webs <- substitute(list(...))
#	#out <- eval(webs)
#	out
#}
#webinput <- webs2array(vazquenc, vazquec)
