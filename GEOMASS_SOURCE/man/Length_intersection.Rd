\name{Length_intersection}
\alias{Length_intersection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate error model of Length intersection
}
\description{
Function to calculate standard deviations of X-Y coordinates determined point using length intersection method. Function calculate X-Y coordinates standard deviations according to standard deviation of measured lengths which is well known for modern total stations. In addition it can also compute X-Y standard deviation of determined point from standard deviations coordinates of points which has been used for calculation (background). If precission of points is not known standard deviations of X-Y coordinates are computed only according to standard deviation of measurement.
}
\usage{
Length_intersection('filepath', choice = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepath}{Adress of the input file. Data in text file are written in it without header, separator is " ", and decimal separator is "." . File consist from two columns and four rows. First and second row is for stand points X-Y coordinates [m] from which the lengths to determined point were measured. Third line is for determined point X-Y coordinates [m] or for lengths [m] measured on stand point in case the user did not compute coordinates by himself. First column for stand point in first line, second for stand point in second line. Four line first column is for the standard deviation of stand points coordinates [mm], if standard deviation is not known type 0 instead(Error model of background will not be computed). Second column is for standard deviation of measured lengths [mm].  }
  \item{choice}{Bool value. FALSE in case user compute coordinates of determined point by himself, if not TRUE has to be choosen. Default value FALSE}
}
\details{

}
\value{
\item{ellipsis}{Data.frame which is returned by function and containg precise characteristics of determined point.}
\item{Output_Length_intersection.csv}{Csv file containing ellipsis data.frame stored in current working directory.}
}
\references{
[1] Skořepa, Zdeněk: Geodézie 4. Publishing ČVUT, Prague 2014,
[2] https://en.wikipedia.org/wiki/Normal_distribution
[3] https://en.wikipedia.org/wiki/Least_squares
}
\author{
%%  ~~who you are~~
}
\note{
Method based on Least square method
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{length_intersection_calculation}},
\code{\link{Elipsa}}
}
\examples{

## Compute error model of determined point

File <- system.file("extdata", "Length_intersection.txt", package = "GEOMASS")
Results <- Length_intersection(File, TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
