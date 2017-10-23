\name{Elipsa}
\alias{Elipsa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute points X-Y coordinates
}
\description{
Function return the X-Y coordinates of the points belonging to ellipse, specified by its parameters.
}
\usage{
Elipsa(a, b, Xp, Yp, om)
}

\arguments{
  \item{a}{ Size of main-axis a}
  \item{b}{ Size of minor-axis b}
  \item{Xp}{ X coordinate of the centre}
  \item{Yp}{ Y coordinate of the centre}
  \item{om}{ Angle between main-axis a and X-axis. Angle must be in gon !}
}
\details{
Note that funtion compute points that belong to ellipse. Function do not plot the ellipse.
}
\value{

Function return matrix of X-Y coordinates
\item{vec[1]}{ X coordinates of the points}
\item{vec[2]}{ Y coordinates of the points}

}
\references{

}
\author{

}
\examples{
## Get ellipse points

points <- Elipsa(1000, 50, 0, 0, 0)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
