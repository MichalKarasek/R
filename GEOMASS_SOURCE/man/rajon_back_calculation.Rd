\name{rajon_back_calculation}
\alias{rajon_back_calculation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute X-Y coordinates using Rajon back method
}
\description{
Compute X-Y coordinates of determined point using Rajon back method.
}
\usage{
rajon_back_calculation(values, gon, choice = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{values}{Text file consist of two columns and three rows. Data in text file are written in it without header, separator is " ", and decimal separator is "." . First two rows of the first column are for X-coordinates [m], and second column for the Y-coordinates [m] of the orientations. Last third line is for measured angle [gon, degree], and distance between determined point and orientation.}
  \item{gon}{ Bool value. True if angles in gon.}
  \item{choice}{Bool value. True to calculate coordinates.}
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
## Compute X-Y coordinates of determined point

File <- system.file("extdata", "Rajon_back_calculation.txt", geopackage = "GEOMASS")
Results <- rajon_back_calculation(File, TRUE, TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
