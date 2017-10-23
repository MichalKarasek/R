\name{rajon_calculation}
\alias{rajon_calculation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute X-Y coordinates of determined points using Rajon method
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rajon_calculation(values, gon, choice)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{values}{Text file consist of two columns and three rows. Data in text file are written in it without header, separator is " ", and decimal separator is "." . First row is for X-Y coordinates [m] of standpoint. Second row for the X-Y coordinates [m] of the orientation. Third line first column is for measured angle [gon] between orientation and determined point, second colums is for length [m] betwwen stand point and determined point.
}
  \item{gon}{Bool value. Angles measured in gons, gon = TRUE(default). Angles measured in degrees, gon = FALSE.
}
  \item{choice}{Bool value. TRUE to calculate coordinates.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{matrix}{Function return matrix cantaining X-Y coordinates of determined point.}
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
# calculate X-Y coordinates of determined point using Rajon method
File <- system.file("extdata", "rajon_calculation.txt", package = "GEOMASS")
Results <- rajon_calculation(File, TRUE, TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
