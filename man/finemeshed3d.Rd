% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/finemeshed3d.R
\name{finemeshed3d}
\alias{finemeshed3d}
\title{finemeshed3d}
\usage{
finemeshed3d(
  x,
  n.breaks1 = 10,
  n.breaks2 = 1000,
  n.breaks3 = 500,
  round.breaks1 = 9,
  round.breaks2 = 9,
  round.breaks3 = 9
)
}
\arguments{
\item{x}{vector or one column matrix to scale. This vector comes from the
vectorization of the 3d array to be transformed. x is named using the
concatenation of the names of the dimension of the array to be transformed,
as shown in the example below.}

\item{n.breaks1}{integer. Number of breaks needed for the first variable
The variable for which classes are in first position when constructing x's
names is the first variable.}

\item{n.breaks2}{integer. Number of breaks needed for the second variable.
The variable for which classes are in second position when constructing x's
names is the second variable.}

\item{n.breaks3}{integer. Number of breaks needed for the third variable.
The variable for which classes are in third position when constructing x's
names is the third variable.}

\item{round.breaks1}{integer. Number of decimals for breaks of the first
variable.}

\item{round.breaks2}{integer. Number of decimals for breaks of the second
variable.}

\item{round.breaks3}{integer. Number of decimals for breaks of the third
variable.}
}
\value{
Returns:
\describe{
\item{finemeshed.array}{Array of dimension n.breaks1 x n.breaks2 x n.breaks3.
The dimension names of finemeshed.array are the breaks created from each
variable and the associated n.breaks. Each value of finemeshed.array is
equal to the value of x indexed by the classes containing the row and column
names of finemeshed.array}
\item{finemeshed.values1}{First variable breaks}
\item{finemeshed.values2}{Second variable breaks}
\item{finemeshed.values3}{Third variable breaks}
}
}
\description{
Function that helps to transform  a vector into a 3 dimensional array (with a fine
mesh). In the implementation of the spiceFP approach, it allows to transform
matrices of coefficients having different dimensions into matrices of the
same dimension in order to perform arithmetic operations. In practice, the
3d array to be transformed is associated with a contingency table, which
implies numerical variables for which classes have been created.
}
\details{
This function is designed to return a 3d fine meshed array and
breaks associated.
In order to obtain a fine mesh, a high number of breaks must be fixed.
}
\examples{
set.seed(4)
count_table<-hist_3d(x = rnorm(1000),
                     y = rnorm( 1000,5,0.1),
                     z = rnorm( 1000,2,1),
                     breaks_x = seq(-4, 4, by =1),
                     breaks_y = seq(2, 8, by =1),
                     breaks_z = seq(-3, 6, by =1))$Hist.Values

df.x<-as.data.frame.table(count_table)
x<-df.x$Freq
names(x)<-paste0(df.x$Var1,"_",df.x$Var2,"_",df.x$Var3)

res.fm3d<- finemeshed3d(x,10,50,100)
dim(res.fm3d$finemeshed.array)
}
