\name{Rajon}
\alias{Rajon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate error model of Rajon.
}
\description{
Function to calculate standard deviations of X-Y coordinates determined point using Rajon method. Function calculate X-Y coordinates standard deviations according to standard deviation of measured length and angle which is well known for modern total stations. In addition it can also compute X-Y standard deviation of determined point from standard deviations coordinates of points which has been used for calculation (background). If precission of points is not known standard deviations of X-Y coordinates are computed only according to standard deviation of measurement.
}
\usage{
Rajon(filepath, choice = FALSE, gon = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepath}{Text file consist of two columns and six rows. Data in text file are written in it without header, separator is " ", and decimal separator is "." . First row is for X-Y coordinates [m] of standpoint. Second row for the X-Y coordinates [m] of the orientation.  Third line is for determined point X-Y coordinates [m] or for measured angle [gon] between orientation and determined point and length [m] between stand point and determined point. Fourth line is for satndard deviation of measured angles [mgon] and measured lengths [mm]. Last fifth and sixth lines are for the standard deviations X-Y coordinates of stand point and orientation.  If precission of points is not known standard deviations of X-Y coordinates are computed only according to standard deviation of measurement (file consist of 4 rows).
}
  \item{choice}{Bool value. FALSE in case user compute coordinates of determined point by himself, if not TRUE has to be choosen. Default value FALSE
}
  \item{gon}{Bool value. If coordinates of determined point need to be calculated set wheter the angle was measured in gons = TRUE(dafault), or degrees = FALSE.
%%     ~~Describe \code{gon} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{ellipsis}{Data.frame which is returned by function and containg precise characteristics of determined point.}
\item{Output_Rajon.csv}{Csv file containing ellipsis data.frame stored in current working directory.}
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
\code{\link{rajon_calculation}},
\code{\link{Elipsa}}
}
\examples{
## Compute error model of Rajon back method

File <- system.file("extdata", "Rajon.txt", package = "GEOMASS")
Result <- Rajon(File, TRUE)
}
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
