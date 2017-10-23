\name{Helmert_transformation}
\alias{Helmert_transformation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Helmert transformation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Function to convert local coordinate systems to standard coordinate Systems within Helmert transformation method and calculate standard deviations of X-Y coordinates determined point using method of least squares.
Function works by comparing the difference in mathematical relations between standard and local coordinate system. Using these comparisons, the core coordinate network is refined based on the calculated transformation key.
Also, the function lists all important statistical values of the counterbalance based on alignment, that is rate of identity and mean absolute difference of coordinate deviation. Function transform new measured points after calculating of transform key.
}
\usage{
Helmert_transformation(filepath)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepath}{
Adress of the input file. Data in text file are written with header, separator is " ", and decimal separator is ".".
The headers of the individual columns are used to recognize the number of identical points and measured points.        X-Y are coordinates defined in the standard system.x-y are coordinates defined in the local system.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{Output_Helmert_transformation.csv}Csv file containing ellipsis data.frame stored in current working directory.
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
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Calculate X-Y coordinates of determined point using length intersection method

File <- system.file("extdata", "helmert.txt", geopackage = "GEOMASS")
Results <- Helmert_transformation(File, TRUE, TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
