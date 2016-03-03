# SUMMARY
# -------
# This file contains a number of useful symbol and function
# definitions that do not fit under any one category. Here is an
# overview of the functions defined in this file:
#
#    caterase(s)
#    sigmoid(x)
#    logit10(x)
#    sigmoid10(x)
#    project.onto.interval(x,a,b)
#    repmat(A,m,n)
#    center.columns(X)
#    mat.sq(A)
#    var1(X)
#    diagsq(X,a)
#    is.binary.factor(x)
#    is.binary(x)
#    binfactor2num(x)
#    binary.from.categorical(x,col.names)
#    none.missing.row(x)
#    computemaf(geno)
#    getgenocounts(geno)
#    check.normal.quantiles(x)
#    roots2(a,b,c)
#
# SYMBOL DEFINITIONS
# ----------------------------------------------------------------------
# Shorthand for machine epsilon.
eps <- .Machine$double.eps

# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
# Output the string using 'cat', then move the cursor back to the
# beginning of the string so that subsequent output will overwrite
# this string.
caterase <- function (s)
  cat(s,rep("\b",nchar(s)),sep="")

# ----------------------------------------------------------------------
# Returns the sigmoid of x. The sigmoid function is also known as
# the logistic link function. It is the inverse of logit(x).
sigmoid <- function (x)
  return(1/(1 + exp(-x)))

# ----------------------------------------------------------------------
# Returns the logit of x (using the base 10 logarithm).
logit10 <- function (x)
  log10((x + eps)/(1 - x + eps))

# ----------------------------------------------------------------------
# Returns the base-10 sigmoid of x. It is the inverse of logit10(x).
sigmoid10 <- function (x)
  return(1/(1 + 10^(-x)))

# ----------------------------------------------------------------------
# Returns x projected onto interval [a,b].
project.onto.interval <- function (x, a, b)
  pmin(b,pmax(a,x))

# ----------------------------------------------------------------------
# Does the same thing as repmat(A,m,n) in MATLAB.
repmat <- function (A,m,n) {
  return(kronecker(matrix(1,m,n),A))
}

# ----------------------------------------------------------------------
# Centers the columns of matrix X so that the entries in each column
# of X add up to zero.
center.columns <- function (X) {
  mu <- matrix(colMeans(X),1,ncol(X))
  X  <- X - repmat(mu,nrow(X),1)
  return(X)
}

# ----------------------------------------------------------------------
# Returns the matrix product X'*X.
mat.sq <- function (A)
  return(t(A) %*% A)

# ----------------------------------------------------------------------
# This is the same as VAR(X,1)' in MATLAB.
var1 <- function (X) {
  n <- nrow(X)
  return(apply(X,2,function(x) (n-1)/n*var(x)))
}

# ----------------------------------------------------------------------
# diagsq(X) returns diag(X'*X).
# diagsq(X,a) returns diag(X'*diag(a)*X).
diagsq <- function (X, a = NULL) {
  if (is.null(a)) {
    n <- nrow(X)
    a <- rep(1,n)
  }

  # Compute y = (X.^2)'*a.
  a <- c(a)
  y <- c(a %*% X^2)
  return(y)
}

# ----------------------------------------------------------------------
# Returns TRUE if x is a factor with exactly 2 levels.
is.binary.factor <- function (x)
  return(is.factor(x) & nlevels(x) == 2)
    
# ----------------------------------------------------------------------
# Returns TRUE if all elements of a numeric vector are either 0 or 1,
# ignoring missing values (NA).
is.binary <- function (x)
  return(is.numeric(x) & sum(!(x == 0 | x == 1),na.rm = TRUE) == 0)

# ----------------------------------------------------------------------
# Convert a factor with exactly 2 levels to a numeric vector with
# values 0 or 1.
binfactor2num <- function (x) {
  if (!is.binary.factor(x))
    stop("Factor must have exactly 2 levels")
  return(as.numeric(x) - 1)
}
  
# ----------------------------------------------------------------------
# Returns a data frame with one column for each level of factor x. The
# columns of the data frame encode the categorical variable with n
# levels as n binary variables.
binary.from.categorical <- function (x, col.names = NULL) {

  # Create a binary factor for each value of the categorical variable.
  d <- list()
  for (value in levels(x))
    d[[value]] <- factor(as.integer(x == value))

  # Convert the list to a data frame, and adjust the column names, if
  # requested.
  d <- data.frame(d,check.names = FALSE)
  if (!is.null(col.names))
    names(d) <- col.names

  # Output the newly created data frame.
  return(d)
}

# ----------------------------------------------------------------------
# For each row of the matrix or data frame, returns true if all the
# entries in the row are provided (not missing).
none.missing.row <- function (x)
    rowSums(is.na(x)) == 0

# ----------------------------------------------------------------------
# Returns the minor allele frequency given a vector of genotypes
# encoded as allele counts.
computemaf <- function (geno) {
  f <- mean(geno,na.rm = TRUE)/2
  return(min(f,1-f))
}

# ----------------------------------------------------------------------
# Return a data frame with three columns, "AA", "AB" and "BB", that
# report the number of homozygous AA, heterozygous and homozygous BB
# genotypes, respectively, at each SNP. The input is a matrix or data
# frame with columns corresponding to SNPs, and rows corresponding to
# samples.
getgenocounts <- function (geno)  {
  counts <- cbind(colSums(geno == 0,na.rm = TRUE),
                  colSums(geno == 1,na.rm = TRUE),
                  colSums(geno == 2,na.rm = TRUE))
  colnames(counts) <- c("AA","AB","BB")
  return(counts)
}

# ----------------------------------------------------------------------
# Check whether the observed quantiles match what we would expect
# under the normal distribution.
check.normal.quantiles <- function (x) {

  # Discard the missing values.
  x <- x[!is.na(x)]
    
  # Transform the observations so that they have zero mean and unit
  # variance.
  x <- (x - mean(x))/sqrt(var(x))
  
  # Create a data frame giving the observed and expected proportion of
  # samples within 1, 2 and 3 standard deviations of the mean.
  return(data.frame(exp = c(pnorm(1) - pnorm(-1),
                            pnorm(2) - pnorm(-2),
                            pnorm(3) - pnorm(-3)),
                    obs = c(mean(-1 < x & x < 1),
                            mean(-2 < x & x < 2),
                            mean(-3 < x & x < 3)),
                    row.names = c("sd1","sd2","sd3")))
}

# ----------------------------------------------------------------------
# ROOTS2(A,B,C) returns solutions X to quadratic A*X^2 + B*X + C = 0.
roots2 <- function (a, b, c) {
  q <- -(b + sign(b)*sqrt(b^2 - 4*a*c))/2
  return(c(q/a,c/q))
}
