# Functions to fit generate two-way structured data

Reciprocal <- function(x){
	n = length(x)
	temp = c()
	for(i in 1:n){
		if(x[i] == 0) temp[i] = 0
		else temp[i] = 1/x[i]
	}
	return(temp)
}

# Compute Laplacian matrix of adjacency matrix A (works for weighted A)
laplacian <- function (A, normalised = F) 
{
    n = dim(A)[1]
    temp = apply(abs(A), 2, sum)
    D = diag(temp, nrow = n)
    temp1 = Reciprocal(sqrt(temp))
    half.D = diag(temp1, nrow = n)
    if (normalised == TRUE) 
        return(half.D %*% (D - A) %*% half.D)
    if (normalised == FALSE) 
        return(D - A)
}

# Generate correlated features
gen.one = function(cs) {
  nclusts = p/cs # cs = clust size
  x = NULL
  for (i in 1:nclusts) {
    tf = rnorm(1)
    g = rnorm(cs - 1, tf*rep(0.35, cs - 1), 1)
    x = c(x, tf, g)
  }
  return(x)
}

