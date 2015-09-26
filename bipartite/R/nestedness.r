
nestedness <-function(m, null.models=TRUE,
                      n.nulls=100, popsize=30, n.ind=7,n.gen=2000,
                      binmatnestout=FALSE
                       )
{

# calculates matrix temperature using the binmatnest programm from Miguel Rodriguez-Girones
# Rodriguez-Girones & Santamaria (2006). A new algorithm to calculate the nestedness
# temperature of presence-absence matrices. Journal of Biogeography 33:924-935.
#
# make sure matrix is a valid one as error proofing in the C++ function does not fully work
# and R may crash...
m<- ifelse(m>0,1,0)   # create a binary matrix
if (popsize<n.ind) n.ind <-popsize- 1 # you cannot pick more individuals then there are in the population...
bmn <- .C("bmn5",
          mat=as.integer(m),   #column dominated....
          n.rows = as.integer(nrow(m)),      # notice swapping of cols
          n.cols = as.integer(ncol(m)),      # and rows
          temperature = as.double(-1.0),
          n.nullmodels = as.integer(n.nulls),
          population.size = as.integer(popsize),
          n.individuals = as.integer(n.ind),
          binmatnestout = as.integer(binmatnestout),
          n.generations = as.integer(n.gen),
          nullmodels = as.integer(null.models),
          p.null1 = as.double(-1.0),
          mean.temp.null1 = as.double(-1.0) ,
          var.temp.null1 = as.double(-1.0),
          p.null2 = as.double(-1.0),
          mean.temp.null2 = as.double(-1.0) ,
          var.temp.null2 = as.double(-1.0),
          p.null3 = as.double(-1.0),
          mean.temp.null3 = as.double(-1.0) ,
          var.temp.null3 = as.double(-1.0),
          pack.order.col = as.integer(rep(-1,ncol(m))),
          pack.order.row = as.integer(rep(-1,nrow(m))),
          PACKAGE="bipartite")
if (min(bmn$pack.order.col,bmn$pack.order.row)>-1) bmn$packed.matrix <- m[bmn$pack.order.row,bmn$pack.order.col]

bmn
}
