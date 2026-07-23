# Generate a data set for the JSM tutorial.
library(tools)
library(susieR)
set.seed(1)
p1      <- 2
min_maf <- 0.03
n_out   <- 100

# Load a version of the GTEx genotype data.
X <- readRDS("../data/Thyroid.FMO2.1Mb.RDS")$X
storage.mode(X) <- "double"
X <- X[,1:1000]

# Filter the SNPs by MAF.
maf <- colMeans(X)/2
maf <- pmin(maf,1 - maf)
j   <- which(maf > min_maf)
X   <- X[,j]

# Dimensions of the genotype matrix we will use:
print(dim(X))

# Select the causal SNPs, and then generate the (sparse) true effects.
p        <- ncol(X)
j        <- sample(p,p1)
b        <- rep(0,p)
names(b) <- colnames(X)
b[j]     <- sample(c(-1,1),p1,replace = TRUE)

# Now simulate y.
n <- nrow(X)
e <- rnorm(n,sd = 0.85)
y <- drop(X %*% b + e)
y <- y/sd(y)

# Compute the "in-sample" LD matrix and an "out-of-sample" LD matrix
# (more accurately, the "out-of-sample" matrix this is an LD matrix
# obtained using a subset of the samples).
R <- cor(X)
i <- sample(n,n_out)
R_out <- cor(X[i,])

# Compute the association statistics (z-scores, p-values).
pfromz <- function (z)
  2*pnorm(-abs(z))
out  <- univariate_regression(X,y)
zhat <- with(out,betahat/sebetahat)
pval <- pfromz(zhat)

# Save the data to an .RData file.
b_true <- b
save(file = "jsm_tutorial_data.RData",
     list = c("y","X","zhat","pval","R","R_out","b_true"))
resaveRdaFiles("jsm_tutorial_data.RData")
