\name{intersection_back_calculation}
\alias{intersection_back_calculation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute X-Y coordinates using Intersection back method
}
\description{
Compute X-Y coordinates of determined point using intersection back method.
}
\usage{
intersection_back_calculation(values, gon, choice = FALSE)
}
\arguments{
  \item{values}{Text file consist of two columns and four rows. Data in text file are written in it without header, separator is " ", and decimal separator is "." . First three rows of the first column are for X-coordinates [m], and second column for the Y-coordinates [m] of the points on which the measurements were done. Last fourth line is for measured angles [gon, degree] between points. First column for angle between point on first line and second line. Second column between point on second line and third line.}
  \item{gon}{ Bool value. True if angles in gon.}
  \item{choice}{Bool value. TRUE to calculate coordinates.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{matrix}{Function return matrix containing X-Y coordinates of determined point.}
}
\references{
[1] http://gis.zcu.cz/studium/gen1/html/ch07s03.html,
[2] http://www.fce.vutbr.cz/ged/stavari/vypocty/protinani.html
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Calculate X-Y coordinates of determined point using length intersection method

File <- system.file("extdata", "Intersection_back_calculation.txt", geopackage = "GEOMASS")
Results <- intersection_back_calculation(File, TRUE, TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
