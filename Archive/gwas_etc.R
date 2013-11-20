#make sure we are sure of the Views methods

#make test grange data

#two from each seqlevel
#carry out Views individually for each one - this is the correct data set
#now try - shuffling the grange object
#what order do our views come out
#remove a seqlevel




library(GWAStools,lib.loc='/g/furlong/Harnett/R/')
library(gdsfmt)
createfn.gds('tmp.gds')


# cteate the GDS file "test.gds"
(f <- createfn.gds("test.gds"))
(n <- add.gdsn(f, "matrix", val=matrix(1:(10*6), nrow=10)))
read.gdsn(index.gdsn(f, "matrix"))
# Apply functions over rows of matrix
tmp <- apply.gdsn(n, margin=1, FUN=function(x) print(x))
tmp <- apply.gdsn(n, margin=1,
                  selection = list(rep(c(TRUE, FALSE), 5), rep(c(TRUE, FALSE), 3)),
                  FUN=function(x) print(x))
# Apply functions over columns of matrix
tmp <- apply.gdsn(n, margin=2, FUN=function(x) print(x))
tmp <- apply.gdsn(n, margin=2,
                  selection = list(rep(c(TRUE, FALSE), 5), rep(c(TRUE, FALSE), 3)),
                  FUN=function(x) print(x))
# close
closefn.gds(f)