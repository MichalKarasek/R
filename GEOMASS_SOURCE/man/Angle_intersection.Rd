\name{Angle_intersection}
\alias{Angle_intersection}
\title{Calculate error model for intersection from Angles}
\description{
Function to calculate standard deviations of X-Y coordinates determined point using intersection from angles. Function calculate X-Y coordinates standard deviations according to standard deviation of measured angles which is well-known for modern total stations. In addition it can also compute X-Y standard deviation of determined point from standard deviations coordinates of points which has been used for calculation (background). If precission of point is not known standard deviations of X-Y coordinates are computed only according to standard deviation of measurement.
}
\usage{
Angle_intersection(filepath, sigmaW, choice = FALSE, gon = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filepath}{Adress of the input file. Data in text file are written in it without header, separator is " ", and decimal separator is "." . File consist from two columns and five rows. First and second column is for stand points X-Y coordinates [m] from which the angles to determined point were measured. Third line is for determined point X-Y coordinates [m] or for angles [gon/degree] measured on stand point in case the user did not compute coordinates by himself. First column for stand point in first line, second for stand point in second line. Fourth line is optional and is for standard deviation of X-Y coordinates point in first line. Fifth line is optional and is for standard deviation of X-Y coordinates point in second line. If standard deviations of determined points are calculated from background both fourth and fifth lines need to be specified.}
  \item{sigmaW}{Standard deviation of angle measurement in mgon.}
  \item{choice}{Bool value. TRUE if coordinates of determined points need to be calculate. FALSE if coordinates were computed by user. Default value set to FALSE}
  \item{gon}{Bool value.If choice set to TRUE, units of data to compute coordinates of determined points need to be specified.TRUE for gons, FALSE for degrees. Default value set to TRUE}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{ellipsis}{Data.frame which is returned by function and containg precise characteristics of determined point.}
\item{Output_Angle_intersection.csv}{Csv file containing ellipsis data.frame stored in current working directory.}
}
\references{
[1] Skořepa, Zdeněk: Geodézie 4. Publishing ČVUT, Prague 2014, [2] https://en.wikipedia.org/wiki/Least_squares
}
\author{
%%  ~~who you are~~
}
\note{
Method based on Least square method.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{angle_intersection_calculation}},
\code{\link{Elipsa}}
}
\examples{
# Angle intersection with computed determined point and standard deviation of angle 0.7 mgon

File <- system.file("extdata", "Angle_intersection.txt", package = "GEOMASS")
Results <- Angle_intersection(File, 0.7, FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
