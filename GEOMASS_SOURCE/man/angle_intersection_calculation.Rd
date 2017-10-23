\name{angle_intersection_calculation}
\alias{angle_intersection_calculation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intersection from angles
}
\description{
Compute X-Y coordinates of determined point using angle intersection method.
}
\usage{
angle_intersection_calculation(values, gon, choice)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{values}{Text file consist of two columns and three rows. Data in text file are written in it without header, separator is " ", and decimal separator is "." . First two rows of the first column are for X-coordinates [m], and second column for the Y-coordinates [m] of the stand points on which the measurements were done. Last third line is for measured angles [gon || degree] between determined point and stand points. First column angle for the stand point with coordinates in first line of textfile, second column angle for second standpoint with coordinates in second line.

}
  \item{gon}{
Bool value. Angles measured in gons, gon = TRUE(default). Angles measured in degrees, gon = FALSE.
}
  \item{choice}{
Bool value. TRUE to calculate coordinates.
}

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
# calculate coordinates X-Y of determined point using angle intersection

File <- system.file("extdata", "Angle_intersection_calculation.txt", package = "GEOMASS")
Results <- angle_intersection_calculation(File, TRUE, TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
