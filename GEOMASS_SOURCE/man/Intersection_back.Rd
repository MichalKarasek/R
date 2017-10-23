\name{Intersection_back}
\alias{Intersection_back}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate error model for back intersection method.
}
\description{
Function to calculate standard deviations of X-Y coordinates determined point using back intersection method. Function calculate X-Y coordinates standard deviations according to standard deviation of measured lengths and angles which are well known for modern total stations. In addition it can also compute X-Y standard deviation of determined point from standard deviations coordinates of points which has been used for calculation (background). If precission of points is not known standard deviations of X-Y coordinates are computed only according to standard deviation of measurement.
}
\usage{
Intersection_back(filepath, choice, gon)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepath}{Adress of the input file. Data in text file are written in it without header, separator is " ", and decimal separator is "."  File consist from two columns and five rows. From first to third line is place  for X-Y coordinates [m] for which the lengths to determined point were measured.
}
  \item{choice}{Bool value. FALSE in case user compute coordinates of determined point by himself, if not TRUE has to be choosen. Default value FALSE}

  \item{gon}{Bool value.If choice set to TRUE, units of data to compute coordinates of determined points need to be specified.TRUE for gons, FALSE for degrees. Default value set to TRUE}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{ellipsis}{Data.frame which is returned by function and containg precise characteristics of determined point.}
\item{Output_Back_intersection.csv}{Csv file containing ellipsis data.frame stored in current working directory.}
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

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{Elipsa}}
  \code{\link{intersection_back_calculation}}
}
\examples{
# Calculate standard deviations of X-Y coordinates determined points using back intersection method.

File <- system.file("extdata", "Intersection_back.txt", package="GEOMASS")
Results <- Intersection_back(File) # Txt file has already calculated X-Y coordinates of determined points by user.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
